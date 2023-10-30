import argparse
import subprocess
from multiprocessing import Pool
from functools import partial
import tempfile

from pprint import pprint

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from glob import glob
from matplotlib import pyplot as plt
import tempfile
from pathlib import Path
import re
import bisect
import subprocess
from collections import namedtuple, defaultdict, Counter
import os
from Bio.Align import substitution_matrices
from itertools import product, combinations
import numpy as np
# from Bio.AlignIO.MafIO import MafWriter, MafIterator
# from Bio.AlignIO.MauveIO import MauveWriter, MauveIterator
from logger import logger
import time
import copy
from tqdm import tqdm
import pyabpoa as pa
import spoa
from intervaltree import IntervalTree
from scipy import stats
import math


FASTA_SUFFIX_LIST = ".fasta, .fas, .fa, .fna, .ffn, .faa, .mpfa, .frn".split(", ")
CHUNK_PREFIX = "chunk"


def validate_xmfa(xmfa_file, parsnp_header=False, gid_to_records=None, gid_to_cid_to_index=None, input_dir):
    """
    Validate records in xmfa file
    """
    print(f"Validating {xmfa_file}")
    seqid_parser = re.compile(r'^cluster(\d+) s(\d+):p(\d+)')
    
    index_to_gid, gid_to_index = parse_xmfa_header(xmfa_file)
    if gid_to_cid_to_index is None or gid_to_records is None:
        gid_to_records, gid_to_cid_to_index = index_input_sequences(xmfa_file, input_dir)
    for lcb in tqdm((AlignIO.parse(xmfa_file , "mauve"))):
        for seq in lcb:
            if parsnp_header:
                aln_len = seq.annotations["end"] - seq.annotations["start"] + 1
                # if " :p" in seq.id:
                #     seq.id = seq.id.replace(" :p", " s1:p")
                #     offset = 1
                #     print(seq.id)
                # else:
                #     offset = 0
                try:
                    cluster_idx, contig_idx, startpos = [int(x) for x in seqid_parser.match(seq.id).groups()]
                except:
                    print(seq.id)

                gid = index_to_gid[seq.name]
                cid = gid_to_cid_to_index[gid][contig_idx]

                if seq.annotations["strand"] == -1:
                    startpos, endpos = startpos - aln_len, startpos
                else:
                    endpos = startpos + aln_len
            else:
                startpos, endpos, strand = seq.annotations["start"], seq.annotations["end"], seq.annotations["strand"]
                aln_len = endpos - startpos
                gid, cid = seq.id.split("#")
            # assert endpos - startpos == len(seq.seq)
            true_str = gid_to_records[gid][cid].seq[startpos:endpos].lower()
            if seq.annotations["strand"] == -1:
                true_str = true_str.reverse_complement()

            true_str = "".join(x if x in "actg" else "n" for x in str(true_str).lower())
            aln_str = "".join(x if x in "actg" else "n" for x in str(seq.seq).lower().replace("-", ""))
            if len(true_str) != len(aln_str) or true_str.lower()[:5] != aln_str[:5]:
                print(startpos, endpos)
                print(seq.name)
                print(seq.id)
                print(seq.annotations)
                print(gid, cid, seq.annotations["strand"])
                print(len(seq.seq.replace("-", "")), aln_len)
                print(true_str[:10], true_str[-10:])
                print(seq.seq.replace("-", "")[:10].lower(), seq.seq.replace("-", "")[-10:].lower())
                print("expected str", true_str)
                print("record str  ", seq.seq.replace("-", "").lower())
                print(seq.seq)
                print()
                print()
                # raise
        

def getIntersection(interval_1, interval_2):
    start = max(interval_1[0], interval_2[0])
    end = min(interval_1[1], interval_2[1])
    if start < end:
        return (start, end)
    return None

def intersect(intervals1, intervals2):
    if len(intervals1) == 0 or len(intervals2) == 0:
        return 
    iter1 = iter(intervals1)
    iter2 = iter(intervals2)

    interval1 = next(iter1)
    interval2 = next(iter2)

    while True:
        intersection = getIntersection(interval1, interval2)
        if intersection:
            yield intersection
            try:
                if intersection[1] == interval1[1]:
                    interval1 = next(iter1)
                else:
                    interval2 = next(iter2)
            except StopIteration:
                return

        try:
            # If first interval starts after the end of the second interval
            # inc the second interval
            while interval1[0] >= interval2[1]:
                interval2 = next(iter2)
             
            # If second interval starts after the end of the first interval
            # inc the first interval
            while interval2[0] >= interval1[1]:
                interval1 = next(iter1)
        except StopIteration:
            return


def get_interval(lcb, seqidx):
    seqid_parser = re.compile(r'^cluster(\d+) s(\d+):p(\d+)/.*')
    for seq in lcb:
        if int(seq.name) == seqidx:
            aln_len = seq.annotations["end"] - seq.annotations["start"] + 1
            cluster_idx, contig_idx, startpos = [int(x) for x in seqid_parser.match(seq.id).groups()]

            if seq.annotations["strand"] == -1:
                startpos, endpos = startpos - aln_len, startpos
            else:
                endpos = startpos + aln_len
                
            return (contig_idx, (startpos, endpos))
        

def cut_overlaps(ilist):
    for i in range(len(ilist) - 1):
        if ilist[i][1] - 1 >= ilist[i+1][0]:
            ilist[i+1] = (ilist[i][1], ilist[i+1][1])

def trim(lcb, ref_cidx_to_intervals, seqidx, cluster_start):
    seqid_parser = re.compile(r'^cluster(\d+) s(\d+):p(\d+)')
    ret_lcbs = []
    empty_seqs = [
        SeqRecord(Seq(""), id=rec.id, name=rec.name, description=rec.description, annotations=copy.deepcopy(rec.annotations))
        for rec in lcb
    ]
    for rec in lcb:
        if int(rec.name) == seqidx:
            aln_len = rec.annotations["end"] - rec.annotations["start"] + 1
            cluster_idx, contig_idx, super_startpos = [int(x) for x in seqid_parser.match(rec.id).groups()]

            if rec.annotations["strand"] == -1:
                super_startpos, super_endpos = super_startpos - aln_len, super_startpos
            else:
                super_endpos = super_startpos + aln_len
                
            try:
                trimmed_intervals = list(intersect(ref_cidx_to_intervals[contig_idx], [(super_startpos, super_endpos)]))
                assert(all(interval in ref_cidx_to_intervals[contig_idx] for interval in trimmed_intervals))
                ref_rec = rec
            except:
                print(trimmed_intervals)
                print((super_startpos, super_endpos))
                print(ref_cidx_to_intervals[contig_idx])
                raise
            break
    else:
        print("Interval not found!")
        print(lcb)
        print(seqidx)
        raise
    
    # ref_psum[i] = number of nucleotides in first i columns
    ref_psum = [0]*(len(ref_rec.seq) + 1)
    for i in range(1, len(ref_rec.seq)+1):
        ref_psum[i] = ref_psum[i-1] + (0 if ref_rec.seq[i-1] == '-' else 1)
        
    # ref_ssum[i] = number of nucleotides in last i columns
    ref_ssum = [0]*(len(ref_rec.seq) + 1)
    for i in range(1, len(ref_rec.seq)+1):
        ref_ssum[i] = ref_ssum[i-1] + (0 if ref_rec.seq[-i] == '-' else 1)
    
    lcb_psum = [[0]*(len(ref_rec.seq) + 1) for _ in range(len(lcb))]
    lcb_ssum = [[0]*(len(ref_rec.seq) + 1) for _ in range(len(lcb))]
    
    for rec_idx, rec in enumerate(lcb):
        psum, ssum = lcb_psum[rec_idx], lcb_ssum[rec_idx]
        # psum[i] = number of nucleotides in first i columns
        for i in range(1, len(ref_rec.seq)+1):
            psum[i] = psum[i-1] + (0 if rec.seq[i-1] == '-' else 1)

        # ssum[i] = number of nucleotides in last i columns
        for i in range(1, len(ref_rec.seq)+1):
            ssum[i] = ssum[i-1] + (0 if rec.seq[-i] == '-' else 1)
    
    for interval_idx, trimmed_interval in (enumerate(trimmed_intervals)):
        lcb_seqs = copy.deepcopy(empty_seqs)
        
        left_bases_trim = trimmed_interval[0] - super_startpos
        right_bases_trim = super_endpos - trimmed_interval[1]
        
        left_cols_to_skip = bisect.bisect_left(ref_psum, left_bases_trim)
        right_cols_to_skip = bisect.bisect_left(ref_ssum, right_bases_trim)
           
        for seq_idx, rec in enumerate(lcb):
            # orig_seq = copy.deepcopy(seq)
            new_rec = lcb_seqs[seq_idx]
            
            aln_len = rec.annotations["end"] - rec.annotations["start"] + 1
            cluster_idx, contig_idx, startpos = [int(x) for x in seqid_parser.match(rec.id).groups()]
            left_bases_trim = 0
            right_bases_trim = 0
                
            left_bases_trim = lcb_psum[seq_idx][left_cols_to_skip]
            right_bases_trim = lcb_ssum[seq_idx][right_cols_to_skip]
    
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
            
        new_lcb = MultipleSeqAlignment(lcb_seqs)
        ret_lcbs.append(new_lcb)
    return ret_lcbs


def copy_header(orig_xmfa, new_xmfa):
    with open(orig_xmfa) as xmfa_in, open(new_xmfa, 'w') as xmfa_out:
        for line in xmfa_in:
            if line[0] == "#":
                xmfa_out.write(line)
            else:
                break


def write_lcb(lcb, out_handle):
    LINESIZE = 80
    for rec in lcb:
        header = f"> {rec.name}:{rec.annotations['start']+1}-{rec.annotations['end']} {'+' if rec.annotations['strand'] == 1 else '-'} {rec.id}\n"
        out_handle.write(header)
        # print(len(rec.seq))
        for i in range(math.ceil(len(rec.seq) / LINESIZE)):
            out_handle.write(str(rec.seq[i*LINESIZE:(i+1)*LINESIZE]) + "\n")
    out_handle.write("=\n")


def combine_header_info(xmfa_list):
    SequenceEntry = namedtuple("SequenceEntry", "index file header length")
    fidx_to_new_idx = {}
    seq_to_idx = {}
    header_to_xmfa = {}
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
                
                entry = SequenceEntry(seqidx, file, header, length)
                # Avoid duplicated headers, i.e. if the reference is duplicated
                if (file, header) not in file_header_pairs:
                    file_header_pairs.add((file, header))
                    fidx_to_new_idx[(xmfa_file, entry.index)] = len(fidx_to_new_idx) + 1
                    seq_to_idx[entry] = fidx_to_new_idx[(xmfa_file, entry.index)]
                    header_to_xmfa[entry.header] = xmfa_file
                line = next(xmfa_in).strip()
            
    return seq_to_idx, fidx_to_new_idx, header_to_xmfa


def write_combined_header(seq_to_idx, cluster_count, xmfa_out_f, header_to_xmfa):
    """
    Writes an XMFA header for the partitioned XMFA files.

    seq_to_idx = {seq_name : idx_in_combined_xmfa}
    cluster_count = number of clusters
    xmfa_out_f = Name of output xmfa
    header_to_xmfa = {xmfa_sequence_header : xmfa_file}
    fidx_to_new_idx 
    """

    #TODO fix duplicated fidx_to_new_idx
    with open(xmfa_out_f, 'w') as xmfa_out:
        xmfa_out.write("#FormatVersion Parsnp v1.1\n")
        xmfa_out.write(f"#SequenceCount {len(seq_to_idx)}\n")
        for idx, entry in enumerate(sorted(seq_to_idx.keys(), key=lambda k: seq_to_idx[k])):
            xmfa_out.write(f"##SequenceIndex {seq_to_idx[entry]}\n")
            xmfa_out.write(f"##SequenceFile {entry.file}\n")
            xmfa_out.write(f"##SequenceHeader {entry.header}\n")
            xmfa_out.write(f"##SequenceLength {entry.length}bp\n")
        
        xmfa_out.write(f"#IntervalCount {cluster_count}\n")

def merge_blocks(aln_xmfa_pairs, fidx_to_new_idx):
    """
    aln_xmfa_pairs = [(aln1, xmfa_file1), ..., (alnn, xmfa_filen)]
    fidx_to_new_idx = {(xmfa_index, contig_idx) : index_in_combined_xmfa}
    """
    # ref is present in every block, so (len(alns[0].seq)-1)*len(alns) + 1 total seqs
    # Set the IDs
    aln, xmfa_file = aln_xmfa_pairs[0]
    combined_seqs = []
    name_to_idx = {}
    xmfa_file_to_col = {}
    for seq in aln[:-1]:
        # This copies the string as well, but we could do it faster by just copying the seq metadat
        new_seq = copy.deepcopy(seq)
        new_seq.name = fidx_to_new_idx[(xmfa_file, int(seq.name))]
        new_seq.id = new_seq.id.split("/")[0]
        xmfa_file_to_col[xmfa_file] = 0
        new_seq.seq = Seq("")
        name_to_idx[new_seq.name] = len(combined_seqs)
        combined_seqs.append(new_seq)
    
    for aln, xmfa_file in aln_xmfa_pairs[1:]:
        for seq in aln[1:-1]:
            new_seq = copy.deepcopy(seq)
            new_seq.name = fidx_to_new_idx[(xmfa_file, int(seq.name))]
            new_seq.id = new_seq.id.split("/")[0]
            xmfa_file_to_col[xmfa_file] = 0
            new_seq.seq = Seq("")
            name_to_idx[new_seq.name] = len(combined_seqs)
            combined_seqs.append(new_seq)
 
    sorted_names = sorted(name_to_idx.keys(), key=lambda x: int(x))
    gap_sequences = defaultdict(list)
    ref_pos = 0
    
    while True:
        in_gap = False
        all_done = True
        for aln, xmfa_file in aln_xmfa_pairs:
            col = xmfa_file_to_col[xmfa_file]
            if col >= len(aln[0].seq) or aln[0].seq[col] == "-":
                in_gap = True
            if col < len(aln[0].seq):
                all_done = False
                
        if all_done:
            break
            
        if not in_gap:
            if len(gap_sequences) != 0:
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

            # Add to alignment
            for aln_idx, (aln, xmfa_file) in enumerate(aln_xmfa_pairs):
                col = xmfa_file_to_col[xmfa_file]
                for rec in aln[0 if aln_idx == 0 else 1:-1]:
                    new_name = fidx_to_new_idx[(xmfa_file, int(rec.name))]
                    new_record = combined_seqs[name_to_idx[new_name]]
                    new_record.seq += rec.seq[col]
                xmfa_file_to_col[xmfa_file] += 1 
        else:
            # We are in a gap for some ref sequence
            # For each alignment, take the slice starting at the current position
            # and ending at the next non-gap reference position for aln_idx, (aln, xmfa_file) in enumerate(aln_xmfa_pairs):
            for aln_idx, (aln, xmfa_file) in enumerate(aln_xmfa_pairs):
                col = xmfa_file_to_col[xmfa_file]
                while col < len(aln[0].seq) and aln[0].seq[col] == "-":
                    for rec in aln[0 if aln_idx == 0 else 1:-1]:
                        if col < len(rec.seq) and rec.seq[col] != "-":
                            new_name = fidx_to_new_idx[(xmfa_file, int(rec.name))]
                            gap_sequences[new_name].append(rec.seq[col])
                    xmfa_file_to_col[xmfa_file] += 1 
                    col = xmfa_file_to_col[xmfa_file]
    
    if len(gap_sequences) != 0:
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
            
        # Add resulting alignment and gaps to lcb
        aln_len = max(len(s) for s in gap_sequences.values())

        for name in sorted_names:
            record = combined_seqs[name_to_idx[name]]
            aligned_seq = "".join(gap_sequences[name]) + "-"*(aln_len - len(gap_sequences[name]))
            record.seq += aligned_seq

        # Clear gap sequences
        gap_sequences = defaultdict(list)
        
    return MultipleSeqAlignment(combined_seqs)


def run_command(cmd):
    """
    Runs the provided command string.
    """
    subprocess.run(cmd, shell=True)
    return

def parse_args():
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

def get_chunked_intervals(partition_output_dir, chunk_labels):
    """
    Returns a mapping with the following structure
    {chunk_id : {contig_idx : [(start1, stop1), ..., (startn, stopn)]}}
    """
    chunk_to_invervaldict = {}
    for chunk_label in chunk_labels:
        orig_xmfa = f"{partition_output_dir}/{CHUNK_PREFIX}-{chunk_label}-out/parsnp.xmfa"
        orig_alns = AlignIO.parse(orig_xmfa, "mauve")
        chunk_to_invervaldict[chunk_label] = defaultdict(list)
        for idx, aln in enumerate((orig_alns)):
            # Get reference idx and the interval of the alignment wrt the reference contig
            ref_cidx, interval = get_interval(aln, 1)
            chunk_to_invervaldict[chunk_label][ref_cidx].append(interval)
            
        # Sort the intervals. (They should be sorted already, but worth double checking)
        # Sometimes parsnp will yield overlapping intervals, so we trim them back
        for ref_cidx, intervals in chunk_to_invervaldict[chunk_label].items():
            intervals.sort()
            cut_overlaps(intervals)

    return chunk_to_invervaldict


def get_intersected_intervals(chunk_to_invervaldict, min_interval_size=10):
    """
    Computes the intersection of all of the intervals in chunk_to_invervaldict
    Each intervaldict is a mapping of (contig_idx : [interval_list])
    """
    first_chunk = list(chunk_to_invervaldict.keys())[0]
    intersected_interval_dict = copy.deepcopy(chunk_to_invervaldict[first_chunk])

    num_clusters = []
    bp_covered = []
    for chunk_label, chunk_intervals in chunk_to_invervaldict.items():
        chunk_aligned_bp = 0
        num_intervals = 0
        for refcidx, intervals in chunk_intervals.items():
            chunk_aligned_bp += sum(b - a for a,b in intervals)
            num_intervals += len(intervals)

        num_clusters.append(num_intervals)
        bp_covered.append(chunk_aligned_bp)
        for refcidx in set(intersected_interval_dict.keys()) | set(chunk_intervals.keys()):
            intersected_interval_dict[refcidx] = list((intersect(intersected_interval_dict[refcidx], chunk_intervals[refcidx])))   
        
    intersected_interval_dict = {
        key: [interval for interval in val if (interval[1] - interval[0]) >= min_interval_size] 
        for key, val in intersected_interval_dict.items()
    }
    intersection_sum = 0
    num_ints = 0
    for refcidx, intervals in intersected_interval_dict.items():
        intersection_sum += sum(b - a for a,b in intervals)
        num_ints += len(intervals)

    print(f"Partition stats: Mean bp covered = {np.mean(bp_covered):.2f}\tMean LCB count = {np.mean(num_clusters):.2f}")
    print(f"After intersection:            {intersection_sum} reference bases over {num_ints} clusters")
    return intersected_interval_dict


def trim_single_xmfa(xmfa_file, intersected_interval_dict):
    chunk_parser = re.compile(f'(.*)-out/parsnp.xmfa')
    chunk_label = chunk_parser.search(xmfa_file).groups()[0]
    orig_alns = AlignIO.parse(xmfa_file, "mauve")
    
    trimmed_xmfa = xmfa_file + ".trimmed"
    copy_header(xmfa_file, trimmed_xmfa)
    
    cluster_start = 1
    with open(trimmed_xmfa, 'a') as trimmed_out:
        for idx, aln in enumerate((orig_alns)):
            new_alns = trim(aln, intersected_interval_dict, seqidx=1, cluster_start=cluster_start)
            cluster_start += len(new_alns)
            for new_aln in new_alns:
                write_lcb(new_aln, trimmed_out)
    
    num_clusters = cluster_start - 1
    return num_clusters


def trim_xmfas(partition_output_dir, chunk_labels, intersected_interval_dict, threads=1):
    """
    Trim all XMFA files so that their LCBs are all wrt the same reference coordinates
    """
    trim_partial = partial(trim_single_xmfa, intersected_interval_dict=intersected_interval_dict)
    orig_xmfa_files = [f"{partition_output_dir}/{CHUNK_PREFIX}-{cl}-out/parsnp.xmfa" for cl in chunk_labels]
    with Pool(threads) as pool:
        num_clusters_per_xmfa = list(tqdm(pool.imap_unordered(trim_partial, orig_xmfa_files), total=len(orig_xmfa_files)))
        #TODO clean up
        if not all(num_clusters_per_xmfa[0] == nc for nc in num_clusters_per_xmfa):
            print("ERROR: One of the partitions has a different number of clusters after trimming...")
            exit(1)
    return num_clusters_per_xmfa[0]


def merge_single_LCB(cluster_idx, aln_xmfa_pairs, tmp_directory, fidx_to_new_idx):
    tmp_xmfa = f"{tmp_directory}/cluster-{cluster_idx}.temp" 
    with open(tmp_xmfa, 'w') as xmfa_out_handle:
        new_aln = merge_blocks(aln_xmfa_pairs, fidx_to_new_idx)
        write_lcb(new_aln, xmfa_out_handle)
    return tmp_xmfa


def merge_single_LCB_star(idx_pairs_tuple, tmp_directory, fidx_to_new_idx):
    merge_single_LCB(idx_pairs_tuple[0], idx_pairs_tuple[1], tmp_directory, fidx_to_new_idx)    


def merge_xmfas(partition_output_dir, chunk_labels, xmfa_out_f, num_clusters, threads=1):
    """
    Take all of the trimmed XMFA files and compile them into a single XMFA
    """

    xmfa_files = [f"{partition_output_dir}/{CHUNK_PREFIX}-{cl}-out/parsnp.xmfa.trimmed"
                  for cl in chunk_labels]
    seq_to_idx, fidx_to_new_idx, header_to_xmfa = combine_header_info(xmfa_files)
    write_combined_header(seq_to_idx, num_clusters, xmfa_out_f, header_to_xmfa)

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
                    total=num_clusters)
            )
    
        with open(xmfa_out_f, 'a') as xmfa_out_handle:
            for cluster_idx in range(num_clusters):
                tmp_xmfa = f"{tmp_directory}/cluster-{cluster_idx}.temp" 
                with open(tmp_xmfa) as tx:
                    xmfa_out_handle.write(tx.read())


if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    query_files = list()
    for arg_path in args.sequences:
        if arg_path.endswith(".txt"):
            with open(arg_path) as input_list_handle:
                for line in input_list_handle:
                    query_files.append(line.strip())
        elif any(arg_path.endswith(suff) for suff in FASTA_SUFFIX_LIST):
            query_files.append(arg_path)
        elif os.path.isdir(arg_path):
            for f in os.listdir(arg_path): 
                if any(f.endswith(suff) for suff in FASTA_SUFFIX_LIST):
                    query_files.append(f"{arg_path}/{f}")

    full_query_list_path = f"{args.output_dir}/input-list.txt"
    with open(full_query_list_path, 'w') as input_list_handle:
        for qf in query_files:
            input_list_handle.write(qf + "\n")

    partition_output_dir = f"{args.output_dir}/partition"
    partition_list_dir = f"{partition_output_dir}/input-lists"
    os.makedirs(partition_list_dir, exist_ok=True)
    run_command(f"split -l {args.partition_size} -a 5 --additional-suffix '.txt' {full_query_list_path} {partition_list_dir}/{CHUNK_PREFIX}-")

    chunk_label_parser = re.compile(f'{CHUNK_PREFIX}-(.*).txt')
    chunk_labels = []
    for partition_list_file in os.listdir(partition_list_dir):
        chunk_labels.append(chunk_label_parser.search(partition_list_file).groups()[0])

    parsnp_commands = []
    for cl in chunk_labels:
        chunk_output_dir = f"{partition_output_dir}/{CHUNK_PREFIX}-{cl}-out"
        os.makedirs(chunk_output_dir, exist_ok=True)
        chunk_command = f"./parsnp -d {partition_list_dir}/{CHUNK_PREFIX}-{cl}.txt -r {args.reference} -o {chunk_output_dir} "
        chunk_command += args.parsnp_flags
        chunk_command += f" > {chunk_output_dir}/parsnp.stdout 2> {chunk_output_dir}/parsnp.stderr"
        parsnp_commands.append(chunk_command)

    print("Running partitions...")
    with Pool(args.threads) as pool:
        list(tqdm(pool.imap_unordered(run_command, parsnp_commands, chunksize=1), total=len(parsnp_commands)))

    print("Computing intersection of all partition LCBs...")
    chunk_to_invervaldict = get_chunked_intervals(partition_output_dir, chunk_labels)
    intersected_interval_dict = get_intersected_intervals(chunk_to_invervaldict)
    
    print("Trimming partitioned XMFAs back to intersected intervals...")
    num_clusters = trim_xmfas(partition_output_dir, chunk_labels, intersected_interval_dict, args.threads)

    os.makedirs(f"{args.output_dir}/merged-out/", exist_ok=True)
    xmfa_out_f =  f"{args.output_dir}/merged-out/parsnp.xmfa"

    print(f"Merging trimmed XMFA files into a single XMFA at {xmfa_out_f}")
    merge_xmfas(partition_output_dir, chunk_labels, xmfa_out_f, num_clusters, args.threads)
