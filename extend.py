from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
from Bio.Align import MultipleSeqAlignment
from glob import glob
import logging
import tempfile
from pathlib import Path
import re
import subprocess
from collections import namedtuple, defaultdict, Counter
import bisect
import os
from Bio.Align import substitution_matrices
from itertools import product, combinations
import numpy as np
from Bio.AlignIO.MafIO import MafWriter, MafIterator
from Bio.AlignIO.MauveIO import MauveWriter, MauveIterator
from logger import logger, TqdmToLogger, MIN_TQDM_INTERVAL
from tqdm import tqdm
import time
import spoa
#%%


# MAX_LEN = 20
MIN_LEN = 10
BASE_LENGTH = 10
length_window = 0.10
window_prop = 0.75

IdPair = namedtuple('IdPair', 'gid cid')
Segment = namedtuple("Segment", "idp start stop strand")

def parse_xmfa_header(xmfa_file):
    index_to_gid = {}
    gid_to_index = {}
    with open(xmfa_file) as parsnp_fd:
        for line in (line.strip() for line in parsnp_fd):
            if line[:2] == "##":
                if line.startswith("##SequenceIndex"):
                    index = line.split(' ')[1]
                elif line.startswith("##SequenceFile"):
                    gid = Path(line.split(' ')[1]).stem
                elif line.startswith("##SequenceLength"):
                    gid_to_index[gid] = index
                    index_to_gid[index] = gid
            elif line[0] != "#":
                break
    return index_to_gid, gid_to_index
                    

def index_input_sequences(xmfa_file, file_list):
    basename_to_path = {}
    for f in file_list:
        basename = str(Path(f).stem)
        basename_to_path[basename] = f
    gid_to_records = {}
    gid_to_cid_to_index = {}
    gid_to_index_to_cid = {}
    with open(xmfa_file) as parsnp_fd:
        for line in (line.strip() for line in parsnp_fd):
            if line[:2] == "##":
                if line.startswith("##SequenceFile"):
                    basename = Path(line.split(' ')[1]).stem
                    p = Path(basename_to_path[basename])
                    gid_to_records[p.stem] = {}
                    gid_to_cid_to_index[p.stem] = {}
                    gid_to_index_to_cid[p.stem] = {}
                    for idx, rec in enumerate(SeqIO.parse(str(p), "fasta")):
                        gid_to_records[p.stem][rec.id] = rec
                        gid_to_cid_to_index[p.stem][rec.id] = idx + 1
                        gid_to_index_to_cid[p.stem][idx + 1] = rec.id
    return gid_to_records, gid_to_cid_to_index, gid_to_index_to_cid


def xmfa_to_covered(xmfa_file, index_to_gid, gid_to_index_to_cid):
    seqid_parser = re.compile(r'^cluster(\d+) s(\d+):p(\d+)/.*')
    idpair_to_segments = defaultdict(list)
    cluster_to_named_segments = defaultdict(list)
    for aln in AlignIO.parse(xmfa_file, "mauve"):
        for seq in aln:
            # Skip reference for now...
            aln_len = seq.annotations["end"] - seq.annotations["start"]
            cluster_idx, contig_idx, startpos = [int(x) for x in seqid_parser.match(seq.id).groups()]
            
            gid = index_to_gid[seq.name]
            
            if seq.annotations["strand"] == -1:
                startpos, endpos = startpos - aln_len, startpos
            else:
                endpos = startpos + aln_len
            
            idp = IdPair(gid, gid_to_index_to_cid[gid][contig_idx])
            seg = Segment(idp, startpos, startpos + aln_len, seq.annotations["strand"])
            idpair_to_segments[idp].append(seg)
            cluster_to_named_segments[cluster_idx].append(seg)
    
    for idp in idpair_to_segments:
        idpair_to_segments[idp] = sorted(idpair_to_segments[idp])
    return idpair_to_segments, cluster_to_named_segments


def run_msa(downstream_segs_to_align, gid_to_records):
    keep_extending = True
    iteration = 0
    seq_lens = [seg.stop - seg.start for seg in downstream_segs_to_align]
    longest_seq = max(seq_lens)
    mean_seq_len = np.mean(seq_lens)
    # if sum(
        # mean_seq_len*(1 - length_window) <= (seg.stop - seg.start) <= mean_seq_len*(1 + length_window) for seg in downstream_segs_to_align) > len(downstream_segs_to_align)*window_prop:
        # base_length = int(mean_seq_len*(1 + length_window))
    # else:
        # base_length = BASE_LENGTH
    
    base_length = BASE_LENGTH
    while keep_extending:
        seqs_to_align = ["A" + (str(
                gid_to_records[seg.idp.gid][seg.idp.cid].seq[seg.start:seg.stop] if seg.strand == 1
                else gid_to_records[seg.idp.gid][seg.idp.cid].seq[seg.start:seg.stop].reverse_complement()
        )[:base_length*(2**iteration)]) for seg in downstream_segs_to_align]
        # print(seqs_to_align)
        consensus, aligned_msa_seqs = spoa.poa(seqs_to_align)
        consensus = consensus.replace("A", "", 1)
        aligned_msa_seqs = [s.replace("A", "", 1) for s in aligned_msa_seqs]
        indel_cols = 0
        for i in range(len(aligned_msa_seqs[0])):
            dash_frac = sum(msa_seq[i] == "-" for msa_seq in aligned_msa_seqs) / len(aligned_msa_seqs)
            if dash_frac > 0.5:
                indel_cols += 1
                if indel_cols >= 10:
                    keep_extending = False
                    if i-indel_cols < MIN_LEN:
                        aligned_msa_seqs = ["" for msa_seq in aligned_msa_seqs]
                    else:
                        aligned_msa_seqs = [msa_seq[:i-indel_cols] for msa_seq in aligned_msa_seqs]
                    break
            else:
                indel_cols = 0
        if BASE_LENGTH*(2**iteration) > longest_seq:
            keep_extending = False
            break
        iteration += 1
    return aligned_msa_seqs


def extend_clusters(xmfa_file, gid_to_index, gid_to_cid_to_index, idpair_to_segments, cluster_to_named_segments, gid_to_records):
    ret_lcbs = []
    seqid_parser = re.compile(r'^cluster(\d+) s(\d+):p(\d+)/.*')
    
    for aln_idx, aln in enumerate(tqdm(
            AlignIO.parse(xmfa_file, "mauve"), 
            total=len(cluster_to_named_segments),
            file=TqdmToLogger(logger, level=logging.INFO),
            mininterval=MIN_TQDM_INTERVAL)):
        # validate_lcb(aln, gid_to_records, parsnp_header=True)
        seq = aln[0]
        cluster_idx, contig_idx, startpos = [int(x) for x in seqid_parser.match(seq.id).groups()]
        segs = cluster_to_named_segments[cluster_idx]
        
        # Get upstream sequences
        assert(segs[0].strand == 1)
        upstream_segs_to_align = []
        downstream_segs_to_align = []
        
        for covered_seg in segs:
            seg_idx = bisect.bisect_left(idpair_to_segments[covered_seg.idp], covered_seg)
            if covered_seg.strand == 1:
                uncovered_start = covered_seg.stop
                if seg_idx == len(idpair_to_segments[covered_seg.idp]) - 1:
                    uncovered_stop = len(gid_to_records[covered_seg.idp.gid][covered_seg.idp.cid].seq)
                else:
                    uncovered_stop = idpair_to_segments[covered_seg.idp][seg_idx+1].start
            else:
                uncovered_stop = covered_seg.start
                if seg_idx == 0:
                    uncovered_start = 0
                else:
                    uncovered_start = idpair_to_segments[covered_seg.idp][seg_idx-1].stop
            downstream_segs_to_align.append(Segment(covered_seg.idp, uncovered_start, uncovered_stop, covered_seg.strand))
           
        aligned_msa_seqs = run_msa(downstream_segs_to_align, gid_to_records)
 
        new_lcb = MultipleSeqAlignment([])
        # Assumes alignments are always in the same order
        new_bp = []
        for seg_idx, (covered_seg, uncovered_seg, aln_str) in enumerate(zip(segs, downstream_segs_to_align, aligned_msa_seqs)):
            # Update segment in idpair_to_segments
            if len(aln_str) < MIN_LEN:
                continue
            new_bp_covered = len(aln_str) - aln_str.count("-")
            # print(f"Extending {covered_seg} by {new_bp_covered}")
            new_bp.append(new_bp_covered)
            new_seq = aln_str
            if covered_seg.strand == 1:
                new_seg = Segment(covered_seg.idp, uncovered_seg.start, uncovered_seg.start + new_bp_covered, covered_seg.strand)
                if new_bp_covered > 0:
                    segs[seg_idx] = Segment(covered_seg.idp, covered_seg.start, new_seg.stop, covered_seg.strand)
            else:
                aln_str = Seq(aln_str).reverse_complement()
                new_seg = Segment(covered_seg.idp, covered_seg.start - new_bp_covered, covered_seg.start, covered_seg.strand)
                if new_bp_covered > 0:
                    segs[seg_idx] = Segment(covered_seg.idp, new_seg.start, covered_seg.stop, covered_seg.strand)
            
            new_record = SeqRecord(
                seq=new_seq,
                id=f"cluster{cluster_idx} s{gid_to_cid_to_index[covered_seg.idp.gid][covered_seg.idp.cid]}:p{new_seg.start if new_seg.strand == 1 else new_seg.stop}",
                name=gid_to_index[covered_seg.idp.gid],
                annotations={"start": new_seg.start, "end": new_seg.stop, "strand": new_seg.strand}
            )

            # if covered_seg.strand == 1:
            new_lcb.append(new_record)

        if len(new_lcb) > 0:
            ret_lcbs.append(new_lcb)
    return ret_lcbs


