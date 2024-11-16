import argparse
import subprocess
import shlex
import re
import bisect
import tempfile
import os
import shutil
import copy
import math
import logging
from pathlib import Path
from multiprocessing import Pool
from functools import partial
from collections import namedtuple, defaultdict, Counter
from typing import List, Tuple, Mapping, TextIO

import numpy as np
import spoa
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from logger import logger, TqdmToLogger, MIN_TQDM_INTERVAL
from tqdm import tqdm


FASTA_SUFFIX_LIST = ".fasta, .fas, .fa, .fna, .ffn, .faa, .mpfa, .frn".split(", ")
CHUNK_PREFIX = "chunk"
XMFASequenceEntry = namedtuple("XMFASequenceEntry", "index file header length")

IntervalType = Tuple[int, int]


def interval_intersection(A: List[IntervalType], B: List[IntervalType]) -> List[IntervalType]:
    """
    Args:
        A:  First list
        B:  Second list
    Returns:
        A list of the intersection of the input lists
    """
    ans = []
    i = j = 0

    while i < len(A) and j < len(B):
        # Let's check if A[i] intersects B[j].
        # lo - the startpoint of the intersection
        # hi - the endpoint of the intersection
        lo = max(A[i][0], B[j][0])
        hi = min(A[i][1], B[j][1])
        if lo < hi:
            ans.append([lo, hi])

        # Remove the interval with the smallest endpoint
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1

    return ans


def get_interval(aln: MultipleSeqAlignment) -> Tuple[int, IntervalType]:
    """
    Get the interval of the first sequence (reference) in a MultipleSeqAlignment object

    Args:
        aln:    A MultipleSeqAlignment object representing an LCB
    Returns:
        A tuple (A, B) where A is the contig idx of the sequence at seqidx, and B is the interval
    """
    seqid_parser = re.compile(r'^cluster(\d+) s(\d+):p(\d+)/.*')
    seq = aln[0]
    aln_len = seq.annotations["end"] - seq.annotations["start"] 
    cluster_idx, contig_idx, startpos = [int(x) for x in seqid_parser.match(seq.id).groups()]

    if seq.annotations["strand"] == -1:
        startpos, endpos = startpos - aln_len, startpos
    else:
        endpos = startpos + aln_len
        
    return (contig_idx, (startpos, endpos))

        

def cut_overlaps(ilist: List[IntervalType]) -> None:
    """
    If two intervals overlap, I1 = (a, b), I2 = (b-10, c)
    then trim the second interval to be (b+1, c)

    Args:
        ilist: List of intervals to be cut
    """
    for i in range(len(ilist) - 1):
        if ilist[i][1] > ilist[i+1][0]:
            ilist[i+1] = (ilist[i][1]+1, ilist[i+1][1])


def trim(aln: MultipleSeqAlignment, 
         ref_cidx_to_intervals: Mapping[int, List[IntervalType]], 
         seqidx: int, 
         cluster_start: int) -> List[MultipleSeqAlignment]:
    """
    Slice the input alignment up such that the coordinates of the output alignments with resepct 
    to sequence seqidx are consistent with intervals of the corresponding contig in the 
    ref_cidx_to_intervals mapping.

    Args:
        aln : MultipleSeqAlignment object to be trimmed
        ref_cidx_to_intervals : {reference_contig_idx : [(start1, end1), ..., (startn, endn)]}
        seqidx : The index of the reference sequence in the Parsnp XMFA header (typically 1)
        cluster_start : The index of the first output cluster in the trimmed alignment
    Returns:
        A list of MultipleSeqAlignment objects
    """
    seqid_parser = re.compile(r'^cluster(\d+) s(\d+):p(\d+)')
    ret_alns = []
    # Store a copy of the SeqRecord objects w/ the correct name, id, annotation dict etc
    empty_seqs = [
        SeqRecord(Seq(""), id=rec.id, name=rec.name, description=rec.description, annotations=copy.deepcopy(rec.annotations))
        for rec in aln
    ]

    # Look for the record in the LCB that represents the reference genome
    for rec in aln:
        if int(rec.name) == seqidx:
            aln_len = rec.annotations["end"] - rec.annotations["start"]
            cluster_idx, contig_idx, super_startpos = [int(x) for x in seqid_parser.match(rec.id).groups()]

            if rec.annotations["strand"] == -1:
                super_startpos, super_endpos = super_startpos - aln_len, super_startpos
            else:
                super_endpos = super_startpos + aln_len
                
            try:
                # super_startpos and super_endpos represent the start and endpoint of this 
                # LCB in the reference. This LCB may be trimmed into multiple LCBs, each of 
                # which is represented by one of the trimmed_intervals
                trimmed_intervals = interval_intersection(
                    ref_cidx_to_intervals[contig_idx], 
                    [(super_startpos, super_endpos)])
                # assert(all(interval in ref_cidx_to_intervals[contig_idx] for interval in trimmed_intervals))
                ref_rec = rec
            except Exception as e:
                logger.critical(e)
                raise e
            break
    else:
        logger.critical("Reference alignment not found!")
        raise
    
    # ref_psum[i] = number of nucleotides in first i columns
    ref_psum = [0]*(len(ref_rec.seq) + 1)
    for i in range(1, len(ref_rec.seq)+1):
        ref_psum[i] = ref_psum[i-1] + (0 if ref_rec.seq[i-1] == '-' else 1)
        
    # ref_ssum[i] = number of nucleotides in last i columns
    ref_ssum = [0]*(len(ref_rec.seq) + 1)
    for i in range(1, len(ref_rec.seq)+1):
        ref_ssum[i] = ref_ssum[i-1] + (0 if ref_rec.seq[-i] == '-' else 1)
    
    # aln_psum[j][i] = number of nucleotides in first i columns of the jth sequence 
    aln_psum = [[0]*(len(ref_rec.seq) + 1) for _ in range(len(aln))]
    # aln_ssum[j][i] = number of nucleotides in last i columns of the jth sequence 
    aln_ssum = [[0]*(len(ref_rec.seq) + 1) for _ in range(len(aln))]
    
    for rec_idx, rec in enumerate(aln):
        psum, ssum = aln_psum[rec_idx], aln_ssum[rec_idx]
        # psum[i] = number of nucleotides in first i columns
        for i in range(1, len(ref_rec.seq)+1):
            psum[i] = psum[i-1] + (0 if rec.seq[i-1] == '-' else 1)

        # ssum[i] = number of nucleotides in last i columns
        for i in range(1, len(ref_rec.seq)+1):
            ssum[i] = ssum[i-1] + (0 if rec.seq[-i] == '-' else 1)
    
    for interval_idx, trimmed_interval in (enumerate(trimmed_intervals)):
        aln_seqs = copy.deepcopy(empty_seqs)
        
        left_bases_trim = trimmed_interval[0] - super_startpos
        right_bases_trim = super_endpos - trimmed_interval[1]
        
        left_cols_to_skip = bisect.bisect_left(ref_psum, left_bases_trim)
        right_cols_to_skip = bisect.bisect_left(ref_ssum, right_bases_trim)
           
        for seq_idx, rec in enumerate(aln):
            # orig_seq = copy.deepcopy(seq)
            new_rec = aln_seqs[seq_idx]
            
            aln_len = rec.annotations["end"] - rec.annotations["start"]
            cluster_idx, contig_idx, startpos = [int(x) for x in seqid_parser.match(rec.id).groups()]
            left_bases_trim = 0
            right_bases_trim = 0
                
            left_bases_trim = aln_psum[seq_idx][left_cols_to_skip]
            right_bases_trim = aln_ssum[seq_idx][right_cols_to_skip]
    
            if rec.annotations["strand"] == -1:
                new_rec.annotations["start"] += right_bases_trim
                new_rec.annotations["end"] -= left_bases_trim
                startpos -= left_bases_trim
            else:
                new_rec.annotations["start"] += left_bases_trim
                new_rec.annotations["end"] -= right_bases_trim
                startpos += left_bases_trim
                
            new_rec.id = f"cluster{cluster_start + interval_idx} s{contig_idx}:p{startpos}"
            if right_cols_to_skip > 0:
                new_rec.seq = rec.seq[left_cols_to_skip:-right_cols_to_skip]
            else:
                new_rec.seq = rec.seq[left_cols_to_skip:]
            
        new_aln = MultipleSeqAlignment(aln_seqs)
        ret_alns.append(new_aln)
    return ret_alns


def copy_header(orig_xmfa: str, new_xmfa: str) -> None:
    """
    Copy header from orig_xmfa to new_xmfa

    Args:
        orig_xmfa: Path to xmfa which has the header to be copied
        new_aln:   Path to new xmfa
    """
    with open(orig_xmfa) as xmfa_in, open(new_xmfa, 'w') as xmfa_out:
        for line in xmfa_in:
            if line[0] == "#":
                xmfa_out.write(line)
            else:
                break


def write_aln_to_fna(aln: MultipleSeqAlignment, out_handle: TextIO) -> None:
    """
    Write the MultipleSeqAlignment to the output handle in fna format.
    
    Args:
        aln :        MultipleSeqAlignment
        out_handle:  An open file handle for the output alignment.
    """
    LINESIZE = 80
    for rec in aln:
        header = f"> {rec.name}:{rec.annotations['start']+1}-{rec.annotations['end']} {'+' if rec.annotations['strand'] == 1 else '-'} {rec.id}\n"
        out_handle.write(header)
        for i in range(math.ceil(len(rec.seq) / LINESIZE)):
            out_handle.write(str(rec.seq[i*LINESIZE:(i+1)*LINESIZE]) + "\n")


def combine_header_info(xmfa_list: List[str]) \
        -> Tuple[Mapping[XMFASequenceEntry, int], Mapping[Tuple[str, int], int]]:
    """
    Combine the headers of multiple partitioned parsnp outputs. Will avoid duplicate headers, i.e.
    if two partitions have a sequence with the same header, only one will be present in the output. 

    Args:
        xmfa_list: A list of xmfa files to be combined.

    Returns:
        A tuple (A, B, C) where
        A: Maps sequence metadata tuple to their index in the combined xmfa file.
        B: Maps (file, original_index) pairs to the index in the combined xmfa file.
    """
    fidx_to_new_idx = {}
    seq_to_idx = {}
    file_header_pairs = set()
    for xmfa_file in xmfa_list:
        with open(xmfa_file) as xmfa_in:
            line = next(xmfa_in).strip()
            while not line.startswith("##"):
                line = next(xmfa_in).strip()
            
            while line.startswith("##"):
                seqidx = int(line.split(" ")[1])
                line = next(xmfa_in).strip()
                file = line.split(" ")[1]
                line = next(xmfa_in).strip()
                header = line.split(" ")[1]
                line = next(xmfa_in).strip()
                length = int(line.split(" ")[1][:-2])
                
                entry = XMFASequenceEntry(seqidx, file, header, length)
                # Avoid duplicated headers, i.e. if the reference is duplicated
                if (file, header) not in file_header_pairs:
                    file_header_pairs.add((file, header))
                    fidx_to_new_idx[(xmfa_file, entry.index)] = len(fidx_to_new_idx) + 1
                    seq_to_idx[entry] = fidx_to_new_idx[(xmfa_file, entry.index)]
                line = next(xmfa_in).strip()
            
    return seq_to_idx, fidx_to_new_idx


def write_combined_header(
    seq_to_idx: Mapping[XMFASequenceEntry, int], 
    cluster_count: int, 
    xmfa_out_f: str) -> None:
    """
    Writes an XMFA header for the merged XMFA of the partitioned XMFA files.

    Args:
        seq_to_idx:         Maps sequence metadata tuple to their index in the combined xmfa file.
        cluster_count:      The number of clusters in the combined output.
        xmfa_out_f:         Name of output xmfa
        fidx_to_new_idx:    A mapping of (xmfa_file, sequence_index) to the index in the combined file. 
    """

    #TODO fix duplicated fidx_to_new_idx
    with open(xmfa_out_f, 'w') as xmfa_out:
        xmfa_out.write("#FormatVersion Parsnp v1.1\n")
        xmfa_out.write(f"#SequenceCount {len(seq_to_idx)}\n")
        for entry in sorted(seq_to_idx.keys(), key=lambda k: seq_to_idx[k]):
            xmfa_out.write(f"##SequenceIndex {seq_to_idx[entry]}\n")
            xmfa_out.write(f"##SequenceFile {entry.file}\n")
            xmfa_out.write(f"##SequenceHeader {entry.header}\n")
            xmfa_out.write(f"##SequenceLength {entry.length}bp\n")
        
        xmfa_out.write(f"#IntervalCount {cluster_count}\n")

def merge_blocks(
    aln_xmfa_pairs: List[Tuple[MultipleSeqAlignment, str]], 
    fidx_to_new_idx: Mapping[Tuple[str, int], int]) -> MultipleSeqAlignment:
    """
    Given a list of MultipleSeqAlignment objects representing the same cluster
    and their originating xmfa files, this function combines them into a single 
    MultipleSeqAlignment object. 

    Merging the alignment requires two important steps:
    (1) The sequence index of the alignments must be updated to reflect their index in the 
        combined XMFA file
    (2) Re-aligning insertion sequences. While there is an equivalence relation for base pairs
        which map to the reference, base pairs which are insertions wrt the reference have 
        no equivalence, and therefore need to be re-aligned to other insertions. This is where 
        we use the SPOA API.

    Args:
        aln_xmfa_pairs:     A list of pairs of MultipleSeqAlignment and the originating xmfa files.
        fidx_to_new_idx:    A mapping of (xmfa_file, sequence_index) to the index in the combined file. 

    Return:
        A MultipleSeqAlignment representing the concatenated alignment of all of the input
        alignments.
    """
    # ref is present in every block, so (len(alns[0].seq)-1)*len(alns) + 1 total seqs
    # Set the IDs
    aln, xmfa_file = aln_xmfa_pairs[0]
    combined_seqs = []
    name_to_idx = {}
    xmfa_file_to_col = {}
    for seq in aln:
        # This copies the string as well, but we could do it faster by just copying the seq metadata
        new_seq = copy.deepcopy(seq)
        new_seq.name = fidx_to_new_idx[(xmfa_file, int(seq.name))]
        new_seq.id = new_seq.id.split("/")[0]
        xmfa_file_to_col[xmfa_file] = 0
        new_seq.seq = Seq("")
        name_to_idx[new_seq.name] = len(combined_seqs)
        combined_seqs.append(new_seq)
    
    for aln, xmfa_file in aln_xmfa_pairs[1:]:
        for seq in aln[1:]:
            new_seq = copy.deepcopy(seq)
            new_seq.name = fidx_to_new_idx[(xmfa_file, int(seq.name))]
            new_seq.id = new_seq.id.split("/")[0]
            xmfa_file_to_col[xmfa_file] = 0
            new_seq.seq = Seq("")
            name_to_idx[new_seq.name] = len(combined_seqs)
            combined_seqs.append(new_seq)
 
    sorted_names = sorted(name_to_idx.keys(), key=lambda x: int(x))
    gap_sequences = defaultdict(list)
    
    all_done = False
    while not all_done:
        in_gap = False
        all_done = True
        for aln, xmfa_file in aln_xmfa_pairs:
            col = xmfa_file_to_col[xmfa_file]
            if col >= len(aln[0].seq) or aln[0].seq[col] == "-":
                in_gap = True
            # If any file has remaining base pairs, we are not done
            if col < len(aln[0].seq):
                all_done = False
                
        if (not in_gap or all_done) and len(gap_sequences) != 0:
            # Perform MSA on gap_sequences
            gap_name_to_idx = {}
            seqs_to_align = []
            for name, bases in gap_sequences.items():
                gap_name_to_idx[name] = len(seqs_to_align)
                seqs_to_align.append("".join(bases))
                
            consensus, aligned_msa_seqs = spoa.poa(seqs_to_align)
            
            # Add resulting alignment and gaps to lcb
            aln_len = max(len(s) for s in aligned_msa_seqs)

            for name in sorted_names:
                record = combined_seqs[name_to_idx[name]]
                if name in gap_name_to_idx:
                    aligned_seq = aligned_msa_seqs[gap_name_to_idx[name]]
                else:
                    aligned_seq = "-"*aln_len 
                record.seq += aligned_seq
                
            # Clear gap sequences
            gap_sequences = defaultdict(list)

        elif not in_gap and not all_done:
            # Add to alignment
            for aln_idx, (aln, xmfa_file) in enumerate(aln_xmfa_pairs):
                col = xmfa_file_to_col[xmfa_file]
                for rec in aln[0 if aln_idx == 0 else 1:]:
                    new_name = fidx_to_new_idx[(xmfa_file, int(rec.name))]
                    new_record = combined_seqs[name_to_idx[new_name]]
                    new_record.seq += rec.seq[col]
                xmfa_file_to_col[xmfa_file] += 1 

        elif not all_done:
            # We are in a gap for some ref sequence
            # For each alignment, take the slice starting at the current position
            # and ending at the next non-gap reference position for aln_idx, (aln, xmfa_file) in enumerate(aln_xmfa_pairs):
            for aln_idx, (aln, xmfa_file) in enumerate(aln_xmfa_pairs):
                col = xmfa_file_to_col[xmfa_file]
                while col < len(aln[0].seq) and aln[0].seq[col] == "-":
                    for rec in aln[0 if aln_idx == 0 else 1:]:
                        if col < len(rec.seq) and rec.seq[col] != "-":
                            new_name = fidx_to_new_idx[(xmfa_file, int(rec.name))]
                            gap_sequences[new_name].append(rec.seq[col])
                    col += 1
                xmfa_file_to_col[xmfa_file] = col

    return MultipleSeqAlignment(combined_seqs)



def run_command(cmd: str, check: bool=True) -> None:
    """
    Runs the provided command string.

    Args:
        cmd:    The command to be run.
        check:  Raise exception if subprocess fails 
    """
    res = subprocess.run(cmd, shell=True, check=check)
    return res.returncode

def parse_args() -> None:
    """
    Parses command line arguments.
    """
    parser = argparse.ArgumentParser(description="""
    Partition Parsnp Parser
    """, formatter_class=argparse.RawTextHelpFormatter)
    #TODO Use lambda to check files and directories
    input_output_args = parser.add_argument_group(title="Input/Output")
    input_output_args.add_argument(
        "-c",
        "--curated",
        action = "store_true",
        help = "(c)urated genome directory, use all genomes in dir and ignore MUMi?")
    input_output_args.add_argument(
        "-d",
        "--sequences",
        type = str,
        nargs = '+',
        required = True,
        help = "A list of files containing genomes/contigs/scaffolds. If the file ends in .txt, each line in the file corresponds to the path to an input file.")
    input_output_args.add_argument(
        "-r",
        "--reference",
        type = str,
        default = "",
        help = "(r)eference genome (set to ! to pick random one from sequence dir)",
        required = True)
    input_output_args.add_argument(
        "-g",
        "--genbank",
        nargs = '+',
        help = "A list of Genbank file(s) (gbk)")
    input_output_args.add_argument(
        "-o",
        "--output-dir",
        type = str,
        default = "parsnp-partition-out")
    input_output_args.add_argument(
        "-n",
        "--partition-size",
        type = int,
        default = 100,
        help = "Maximum size of individual partition")
    input_output_args.add_argument(
        "-t",
        "--threads",
        type = int,
        default = 1,
        help = "Maximum number of partitioned parsnp runs to execute in parallel")
    input_output_args.add_argument(
        "--parsnp-flags",
        type = str,
        default = "",
        help = "Flags other than -d, -r, -g,  and -o to pass to each of the partitioned parsnp runs")

    return parser.parse_args()
####################################################################################################

def get_chunked_intervals(partition_output_dir: str, chunk_labels: List[str])\
        -> Mapping[str, Mapping[int, List[IntervalType]]]:
    """
    Returns the aligned reference intervals for each partition.

    Args:
        partition_output_dir:   The output directory of the partitioned Parsnp runs.
        chunk_labels:           A list of the partition labels/ids.

    Returns:
        A mapping which maps partition ids to the reference intervals for that partition.
            {chunk_label: {contig_idx: [(s1, e1), ... (sn, en)]}}
    """
    chunk_to_invervaldict = {}
    for chunk_label in chunk_labels:
        orig_xmfa = f"{partition_output_dir}/{CHUNK_PREFIX}-{chunk_label}-out/parsnp.xmfa"
        orig_alns = AlignIO.parse(orig_xmfa, "mauve")
        chunk_to_invervaldict[chunk_label] = defaultdict(list)
        for aln in orig_alns:
            # Get reference idx and the interval of the alignment wrt the reference contig
            ref_cidx, interval = get_interval(aln)
            chunk_to_invervaldict[chunk_label][ref_cidx].append(interval)
            
        # Sort the intervals. (They should be sorted already, but worth double checking)
        # Sometimes parsnp will yield overlapping intervals, so we trim them back
        for ref_cidx, intervals in chunk_to_invervaldict[chunk_label].items():
            intervals.sort()
            cut_overlaps(intervals)

    return chunk_to_invervaldict


def get_intersected_intervals(
    chunk_to_invervaldict: Mapping[str, Mapping[int, List[IntervalType]]],
    min_interval_size: int=10)\
        -> Mapping[str, Mapping[int, List[IntervalType]]]:
    """
    Returns the intersection of all of the intervals in chunk_to_invervaldict.
    Args:
        chunk_to_invervaldict:  A mapping of partition IDs to mappings of contigs to intervals
        min_interval_size:      Minimum interval size to retain. Smaller intervals will be dropped.

    Returns:
        A mapping which maps reference contigs to the intersected intervals for that contig.
    """
    first_chunk = list(chunk_to_invervaldict.keys())[0]
    intersected_interval_dict = copy.deepcopy(chunk_to_invervaldict[first_chunk])

    num_clusters = []
    bp_covered = []
    for chunk_intervals in chunk_to_invervaldict.values():
        chunk_aligned_bp = 0
        num_intervals = 0
        for refcidx, intervals in chunk_intervals.items():
            chunk_aligned_bp += sum(b - a for a,b in intervals)
            num_intervals += len(intervals)

        num_clusters.append(num_intervals)
        bp_covered.append(chunk_aligned_bp)
        for refcidx in set(intersected_interval_dict.keys()) | set(chunk_intervals.keys()):
            intersected_interval_dict[refcidx] = interval_intersection(
                intersected_interval_dict[refcidx], 
                chunk_intervals[refcidx])
        
    intersected_interval_dict = {
        key: [interval for interval in val if (interval[1] - interval[0]) >= min_interval_size] 
        for key, val in intersected_interval_dict.items()
    }
    intersection_sum = 0
    num_ints = 0
    for refcidx, intervals in intersected_interval_dict.items():
        intersection_sum += sum(b - a for a,b in intervals)
        num_ints += len(intervals)

    logger.info(f"Partition stats: Mean bp covered = {np.mean(bp_covered):.2f}\tMean LCB count = {np.mean(num_clusters):.2f}")
    logger.info(f"After intersection:            {intersection_sum} reference bases over {num_ints} clusters")
    return intersected_interval_dict


def trim_single_xmfa(
    xmfa_file: str, 
    intersected_interval_dict: Mapping[int, List[IntervalType]]) -> int:
    """
    Given an input xmfa file, creates a new xmfa file with the ".trimmed" extension such that
    f"{xmfa_file}.trimmed is an xmfa file with alignments that correspond to the intervals in
    the intersected_interval_dict.

    Args:
        xmfa_file:                  XMFA file to be trimmed.
        intersected_interval_dict:  A map of reference contig indicies to the intervals. 

    Returns:
        The number of clusters in the trimmed XMFA file.
    """
    orig_alns = AlignIO.parse(xmfa_file, "mauve")
    
    trimmed_xmfa = xmfa_file + ".trimmed"
    copy_header(xmfa_file, trimmed_xmfa)
    
    cluster_start = 1
    with open(trimmed_xmfa, 'a') as trimmed_out:
        for aln in orig_alns:
            new_alns = trim(aln, intersected_interval_dict, seqidx=1, cluster_start=cluster_start)
            cluster_start += len(new_alns)
            for new_aln in new_alns:
                write_aln_to_fna(new_aln, trimmed_out)
                trimmed_out.write("=\n")
    
    num_clusters = cluster_start - 1
    return (xmfa_file, num_clusters)


def trim_xmfas(
    partition_output_dir: str, 
    chunk_labels: List[str], 
    intersected_interval_dict: Mapping[int, List[IntervalType]],
    threads: int=1) -> int:
    """
    Trim all XMFA files so that their LCBs are all wrt the same reference coordinates
    
    Args:
        partition_output_dir:       The directory containing the partitioned outputs.
        chunk_labels:               The partition ids.
        intersected_interval_dict:  A map of reference contig indicies to the intervals. 
        threads:                    The number of threads to use.

    Returns:
        The total number of clusters in each XMFA.
    """
    trim_partial = partial(trim_single_xmfa, intersected_interval_dict=intersected_interval_dict)
    orig_xmfa_files = [f"{partition_output_dir}/{CHUNK_PREFIX}-{cl}-out/parsnp.xmfa" for cl in chunk_labels]
    with Pool(threads) as pool:
        num_clusters_per_xmfa = list(tqdm(
            pool.imap_unordered(trim_partial, orig_xmfa_files), 
            total=len(orig_xmfa_files), 
            file=TqdmToLogger(logger,level=logging.INFO),
            mininterval=MIN_TQDM_INTERVAL))
        #TODO clean up
        if not all(num_clusters_per_xmfa[0][1] == nc for xmfa, nc in num_clusters_per_xmfa):
            logger.critical("One of the partitions has a different number of clusters after trimming...")
            raise
    return num_clusters_per_xmfa[0][1]


def merge_single_LCB(
    cluster_idx: int, 
    aln_xmfa_pairs: List[Tuple[MultipleSeqAlignment, str]], 
    tmp_directory: str, 
    fidx_to_new_idx: Mapping[Tuple[str, int], int]) -> str:
    """
    Merges the alignments in aln_xmfa_pairs into a single XMFA alignment entry. It creates an
    output file in the provided directory representing this alignment.

    Args:
        cluster_idx:        The cluster index of the input alignments.
        aln_xmfa_pairs:     The input alingments and their xmfa files.
        tmp_directory:      A temporary directory to store output in.
        fidx_to_new_idx:    A mapping of (xmfa_file, seqidx) tuples to the idx in the output alignment

    Returns:
        The path to the output alignment file.
    """
    tmp_xmfa = f"{tmp_directory}/cluster-{cluster_idx}.temp" 
    with open(tmp_xmfa, 'w') as xmfa_out_handle:
        new_aln = merge_blocks(aln_xmfa_pairs, fidx_to_new_idx)
        write_aln_to_fna(new_aln, xmfa_out_handle)
    return tmp_xmfa


def merge_single_LCB_star(idx_pairs_tuple, tmp_directory, fidx_to_new_idx):
    """
    A helper function for merging LCBs. See the merge_single_LCB documentation.
    """
    merge_single_LCB(idx_pairs_tuple[0], idx_pairs_tuple[1], tmp_directory, fidx_to_new_idx)    


def merge_xmfas(
    output_dir: str,
    partition_output_dir: str, 
    chunk_labels: List[str], 
    num_clusters: int, 
    threads: int=1,
    write_blocks: bool=False) -> None:
    """
    Take all of the trimmed XMFA files and compile them into a single XMFA

    Args:
        partition_output_dir: Directory of partitioned parsnp runs.
        chunk_labels:           IDs of partitions.
        xmfa_out_f:             Path to merged output XMFA file.
        num_clusters:           Number of clusters in merged XMFA file.
        threads:                Number of threads to use.
    """

    xmfa_out_f = f"{output_dir}/parsnp.xmfa" 
    xmfa_files = [f"{partition_output_dir}/{CHUNK_PREFIX}-{cl}-out/parsnp.xmfa.trimmed"
                  for cl in chunk_labels]
    seq_to_idx, fidx_to_new_idx = combine_header_info(xmfa_files)
    write_combined_header(seq_to_idx, num_clusters, xmfa_out_f)

    alns = [AlignIO.parse(f, "mauve") for f in xmfa_files]

    with tempfile.TemporaryDirectory() as tmp_directory:
        merge_single_LCB_star_partial = partial(
                merge_single_LCB_star,
                tmp_directory=tmp_directory,
                fidx_to_new_idx=fidx_to_new_idx)

        pairs_list = (
            [(next(alns[j]), xmfa_files[j]) for j in range(len(xmfa_files))]
            for i in range(num_clusters)
        )
        
        with Pool(threads) as pool:
            tmp_xmfas = list(tqdm(
                    pool.imap_unordered(merge_single_LCB_star_partial, enumerate(pairs_list)), 
                    total=num_clusters,
                    file=TqdmToLogger(logger,level=logging.INFO),
                    mininterval=MIN_TQDM_INTERVAL)
            )
    
        with open(xmfa_out_f, 'a') as xmfa_out_handle:
            for cluster_idx in range(num_clusters):
                tmp_xmfa = f"{tmp_directory}/cluster-{cluster_idx}.temp" 
                with open(tmp_xmfa) as tx:
                    xmfa_out_handle.write(tx.read())
                if write_blocks:
                    Path(f"{output_dir}/blocks/b{cluster_idx}").mkdir(parents=True)
                    shutil.move(tmp_xmfa, f"{output_dir}/blocks/b{cluster_idx}/seq.fna")


