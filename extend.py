from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
from Bio.Align import MultipleSeqAlignment
from glob import glob
import tempfile
from pathlib import Path
import re
import subprocess
from collections import namedtuple, defaultdict, Counter
import os
from Bio.Align import substitution_matrices
from itertools import product, combinations
import numpy as np
from Bio.AlignIO.MafIO import MafWriter, MafIterator
from Bio.AlignIO.MauveIO import MauveWriter, MauveIterator
from logger import logger
import pyabpoa as pa
import scipy
from matplotlib import pyplot as plt
import time
from tqdm import tqdm
from pprint import pprint
#%%

def get_ani_cutoff(sequences, cutoff=0):
    n = max(len(s) for s in sequences)
    matches = 0
    max_matches = 0
    for i in range(n):
        pileup = Counter(s[i] for s in sequences)
        nongaps = sum(pileup[c] for c in pileup if c not in "-")
        matches += sum(pileup[c]*(pileup[c]-1) for c in pileup if c in "ATGCatgc")
        matches += sum(pileup[ambig] * (sum(pileup[c] for c in pileup if c in "ATGCatgcNn") - 1) for ambig in "Nn")
        max_matches += nongaps * (len(sequences) - 1)
        if matches/max_matches < cutoff:
            return (matches, max_matches), i
    return (matches, max_matches), n


#%%
def get_sequence_data(reference, sequence_files, index_files=False):
    fname_to_seqrecord = {}
    fname_contigid_to_length = {}
    fname_contigidx_to_header = {}
    idx = 0
    for fasta_f in [reference] + sequence_files:
        fname = Path(fasta_f).name
        if index_files:
            fname_to_seqrecord[fname] = SeqIO.index(fasta_f, 'fasta')
            for contig_id, record in fname_to_seqrecord[fname].items():
                fname_contigid_to_length[fname, contig_id] = len(record)
        for ctg_idx, record in enumerate(SeqIO.parse(fasta_f, 'fasta')):
            idx += 1
            record.description = record.description.replace(" ", "^")
            fname_contigidx_to_header[fname, ctg_idx] = record.description
            fname_contigid_to_length[fname, record.description] = len(record)
    return fname_contigid_to_length, fname_contigidx_to_header, fname_to_seqrecord

#%%
# Writes the parsnp output to a MAF file containing. While there are methods
# in HarvestTools that do this, they don't account for multiple contigs,
# which is a future addition we would like to have.
def xmfa_to_maf(xmfa_file, output_maf_file, fname_contigidx_to_header, fname_contigid_to_length):
    index_to_fname = {}
    fname_header_to_gcontigidx = {}
    # Avoid loading the whole alignment in memory
    with open(xmfa_file) as parsnp_fd, tempfile.SpooledTemporaryFile(mode='w') as fasta_out:
        for line in (line.strip() for line in parsnp_fd):
            if line[:2] == "##":
                if line.startswith("##SequenceIndex"):
                    index = int(line.split(' ')[1])
                elif line.startswith("##SequenceFile"):
                    filename = Path(line.split(' ')[1]).name
                    index_to_fname[index] = filename
                elif line.startswith("##SequenceHeader"):
                    header = line.split('>')[1].replace(" ", "^")
                    fname_header_to_gcontigidx[filename, header] = index
            elif line[0] != '=':
                fasta_out.write(line + "\n")
        fasta_out.seek(0)
        xmfa_seqs = SeqIO.parse(fasta_out, "fasta")


        curr_cluster_seqs = []
        curr_cluster_idx = 0
        with open(output_maf_file, 'w') as maf_out:
            writer = MafWriter(maf_out)
            writer.write_header()
        header_parser = re.compile(r"^\s?(\d+):(\d+)-(\d+)\s+([+,-])\s+cluster(\d+)\s+s?(\d+)?:p(\d+)")

        # Stores coordinates of aligned regions on the FORWARD strand
        fname_to_contigid_to_coords = defaultdict(lambda: defaultdict(list))
        for seqrecord in xmfa_seqs:
            full_header = seqrecord.description
            # print(full_header)
            groups = header_parser.search(full_header).groups()
            file_idx, gstart, gend, strand, cluster_idx, contig_idx, position_in_contig = groups
            cluster_idx = int(cluster_idx)
            if cluster_idx != curr_cluster_idx:
                if len(curr_cluster_seqs) != 0:
                    with open(output_maf_file, 'a') as maf_out:
                        writer = MafWriter(maf_out)
                        msa_record = MultipleSeqAlignment(curr_cluster_seqs)
                        msa_record._annotations={"pass": curr_cluster_idx}
                        writer.write_alignment(msa_record)
                    curr_cluster_seqs = []
            curr_cluster_idx = cluster_idx

            # There are some Parsnp indexing situations that need to be corrected.
            if not contig_idx:
                contig_idx = 0
                start_coord = int(position_in_contig) - 1
            else:
                start_coord = int(position_in_contig)
                contig_idx = int(contig_idx) - 1
            if int(file_idx) == 1 and contig_idx != 0:
                start_coord -= 1

            # Add annotations to the seqrecord in order to pass to mafwriter
            fname = index_to_fname[int(file_idx)]
            contig_id = fname_contigidx_to_header[fname, contig_idx]
            # print(contig_id)
            src_size = fname_contigid_to_length[fname, contig_id]
            size = len(str(seqrecord.seq).replace("-", ''))
            seqrecord.id = f"{fname}:<-file:contig->:{contig_id}"
            seqrecord.annotations["strand"] = 1 if strand == "+" else -1
            seqrecord.annotations["srcSize"] = src_size
            seqrecord.annotations["size"] = size
            # Parsnp coordinates are always on the fwd strand, even for revcomp hits
            # we adjust for that here
            seqrecord.annotations["start"] = start_coord if strand == "+" else src_size - start_coord
            curr_cluster_seqs.append(seqrecord)

            # In order to keep track of inter-cluster regions, we also would like
            # to have the coordinates of all aligned clusters on the fwd strand
            if strand == "-":
                end_coord = start_coord
                start_coord = start_coord - size
            else:
                end_coord = start_coord + size
            fname_to_contigid_to_coords[fname][contig_id].append((start_coord, end_coord, strand, cluster_idx))

        # Write remaining sequences
        if len(curr_cluster_seqs) != 0:
            with open(output_maf_file, 'a') as maf_out:
                writer = MafWriter(maf_out)
                msa_record = MultipleSeqAlignment(curr_cluster_seqs)
                msa_record._annotations={"pass": curr_cluster_idx}
                writer.write_alignment(msa_record)
    return fname_to_contigid_to_coords, fname_header_to_gcontigidx
#%%
def write_intercluster_regions(input_sequences, cluster_directory, fname_to_contigid_to_coords):
    clusterdir_to_adjacent_clusters = defaultdict(set)
    fname_contigid_to_cluster_dir_to_length = defaultdict(lambda: defaultdict(int))
    fname_contigid_to_cluster_dir_to_adjacent_cluster = defaultdict(dict)
    for f in glob(f"{cluster_directory}/*.fa*"):
        os.remove(f)
    for fasta_f in input_sequences:
        records = SeqIO.parse(fasta_f, 'fasta')
        fname = Path(fasta_f).name
        for record in records:
            record.description = record.description.replace(" ", "^")
            coords = [(0, 0, "+", "START CAP")] + sorted(fname_to_contigid_to_coords[fname][record.description]) + [(len(record), len(record), "+", "END CAP")]
            for idx, (start, end, strand, cluster_idx) in enumerate(coords[1:-1]):
                idx += 1
                prev_end = coords[idx-1][1]
                next_start = coords[idx+1][0]
                fwd, rev = "right", "left"
                seq = record.seq
                if strand == "-":
                    fwd, rev = rev, fwd
                    # seq = seq.reverse_complement()
                fname_contigid_to_cluster_dir_to_length[(fname, record.description)][(cluster_idx, rev)] = start - prev_end
                fname_contigid_to_cluster_dir_to_length[(fname, record.description)][(cluster_idx, fwd)] = next_start - end
                fname_contigid_to_cluster_dir_to_adjacent_cluster[(fname, record.description)][(cluster_idx, rev)] = (coords[idx-1][3], "right" if coords[idx-1][2] == "+" else "left")
                fname_contigid_to_cluster_dir_to_adjacent_cluster[(fname, record.description)][(cluster_idx, fwd)] = (coords[idx+1][3], "left" if coords[idx+1][2] == "+" else "right")
                clusterdir_to_adjacent_clusters[(cluster_idx, fwd)].add((coords[idx+1][3], "left" if coords[idx+1][2] == "+" else "right"))
                clusterdir_to_adjacent_clusters[(cluster_idx, rev)].add((coords[idx-1][3], "right" if coords[idx-1][2] == "+" else "left"))
                seq_to_write = SeqRecord(seq[prev_end:start].reverse_complement(), id=f"{fname}:<-file:contig->:{record.description}", description="")
                with open(f"{cluster_directory}/cluster{cluster_idx}-{rev}.fasta", 'a') as out_f:
                    SeqIO.write([seq_to_write], out_f, "fasta")
                seq_to_write = SeqRecord(seq[end:next_start], id=f"{fname}:<-file:contig->:{record.description}", description="")
                with open(f"{cluster_directory}/cluster{cluster_idx}-{fwd}.fasta", 'a') as out_f:
                    SeqIO.write([seq_to_write], out_f, "fasta")
    return fname_contigid_to_cluster_dir_to_length, fname_contigid_to_cluster_dir_to_adjacent_cluster
#%%
def get_new_extensions(cluster_files, match_score, mismatch_score, gap_penalty):
    fname_parser = re.compile(r"^.*/cluster(\d+)-(.*).fasta")
    clusterdir_len = defaultdict(lambda: defaultdict(int))
    clusterdir_expand = {}
    for cluster_f in cluster_files:
        cluster_idx, direction = fname_parser.match(cluster_f).groups()
        cluster_idx = int(cluster_idx)
        sequences = list(SeqIO.parse(cluster_f, 'fasta'))
        lengths = Counter(len(s) for s in sequences)
        maxlen = max(lengths)
        align_all = True
        if maxlen != 0:
            bins = np.histogram([len(s) for s in sequences], bins=20)[0]
            if max(bins)/len(sequences) < 0.99 or lengths.most_common(1)[0][0] == 0:
                align_all = False

        if not align_all:
            idx = 0
            score_list = [0]
            clusterdir_len[(cluster_idx, direction)] = max(lengths)
            while idx < max(lengths):
                bases = Counter(seq[idx].upper() if len(seq) > idx else "-" for seq in sequences if seq)
                score = 0
                for n1, n2 in combinations(bases.keys(), r=2):
                    if n1 == "-" or n2 == "-":
                        score += gap_penalty * bases[n1] * bases[n2]
                    else:
                        score += mismatch_score * bases[n1] * bases[n2]
                for nuc in bases.keys():
                    if nuc == "-":
                        continue
                    else:
                        score += match_score * bases[nuc] * (bases[nuc] - 1) / 2
                score_list.append(score_list[-1] + score)
                idx += 1

            extend_by = np.argmax(score_list)
        else:
            extend_by = lengths.most_common(1)[0][0]
        clusterdir_expand[(cluster_idx, direction)] = extend_by
    return clusterdir_expand, clusterdir_len

#%%

def write_xmfa_header(extended_xmfa_file, fname_header_to_gcontigidx, fname_contigid_to_length, num_clusters):

    with open(extended_xmfa_file, 'w') as xmfa_out:#, open(extended_xmfa_file + ".maf", 'w') as maf_out:
        # maf_writer = MafWriter(maf_out)
        # maf_writer.write_header()
        ### Write header
        xmfa_out.write("#FormatVersion Parsnp v1.1\n")
        xmfa_out.write(f"#SequenceCount {len(fname_header_to_gcontigidx)}\n")
        for (fname, header), global_idx in fname_header_to_gcontigidx.items():
            xmfa_out.write(f"##SequenceIndex {global_idx}\n")
            xmfa_out.write(f"##SequenceFile {fname}\n")
            xmfa_out.write(f"##SequenceHeader >{header.replace('^', ' ')}\n")
            xmfa_out.write(f"##SequenceLength {fname_contigid_to_length[(fname, header)]}\n")
        xmfa_out.write(f"##IntervalCount {num_clusters}\n")


def write_xmfa_cluster(extended_xmfa_file, msa_records, fname_header_to_gcontigidx):
    with open(extended_xmfa_file, 'a') as xmfa_out:#, open(extended_xmfa_file + ".maf", 'a') as maf_out:
        # maf_writer = MafWriter(maf_out)
        ### Write MSA record
        header_parser = re.compile(r"(.+):<-file:contig->:(.+)")
        for msa_record in msa_records:
            cluster = msa_record.annotations['cluster']
            # maf_writer.write_alignment(msa_record)
            for seq_record in msa_record:
                fname, contig_id = header_parser.match(seq_record.id).groups()

                strand = "-" if seq_record.annotations["strand"] == -1 else "+"
                start = seq_record.annotations['start']
                size = seq_record.annotations["size"]
                if strand == "+" or True:
                    seq_record.id = f"{fname_header_to_gcontigidx[(fname, contig_id)]}:{start}-{start + size - 1}"
                else:
                    seq_record.id = f"{fname_header_to_gcontigidx[(fname, contig_id)]}:{start - size}-{start}"
                seq_record.description = f"{strand} cluster{cluster} s1:p{start}"
            SeqIO.write(msa_record, xmfa_out, "fasta")
        xmfa_out.write("=\n")


# fname_to_contigid_to_extended = defaultdict(lambda: defaultdict(list))
# fname_to_contigid_to_coords[fname][contig_id].append((seqrecord.annotations["start"], seqrecord.annotations["start"] + seqrecord.annotations["size"], strand, cluster_idx))
def write_extended_xmfa(
        original_maf_file,
        extended_maf_file,
        cluster_directory,
        clusterdir_expand,
        clusterdir_len,
        fname_contigid_to_cluster_dir_to_length,
        fname_contigid_to_cluster_dir_to_adjacent_cluster,
        fname_header_to_gcontigidx,
        fname_contigid_to_length,
        cutoff,
        indel_cutoff,
        cpu_count=1): # cutoff is the min ANI 
    header_parser = re.compile(r"(.+):<-file:contig->:(.+)")
    clusters = set()
    for fname, contigid in fname_contigid_to_cluster_dir_to_length.keys():
        for cluster, _ in fname_contigid_to_cluster_dir_to_length[(fname, contigid)].keys():
            clusters.add(cluster)
    write_xmfa_header(extended_maf_file, fname_header_to_gcontigidx, fname_contigid_to_length, len(clusters))
    with open(original_maf_file, 'r') as maf_in:
        maf_iterator = MafIterator(maf_in)
        # maf_writer = MafWriter(maf_out)
        # maf_writer = MauveWriter(maf_out)
        # maf_writer.write_header()
        old_matches = 0
        old_max_matches = 0
        new_matches = 0
        new_max_matches = 0
        old_nucs_aligned = 0
        new_nucs_aligned = 0
        for idx, msa_record in (enumerate(maf_iterator)):
            fname, contig_id = header_parser.match(msa_record[0].id).groups()
            cluster_idx = int(msa_record._annotations["pass"])
            (record_matches, record_maxmatches), _ = get_match_maxmatch([record.seq for record in msa_record])
            old_matches, old_max_matches = record_matches + old_matches, record_maxmatches + old_max_matches
            old_nucs_aligned += sum(len(record.seq) for record in msa_record)
            for direction in ("right", "left"):
                expand_by = clusterdir_expand[(cluster_idx, direction)]
                flanks = []
                inter_cluster_sequences = SeqIO.index(f"{cluster_directory}/cluster{cluster_idx}-{direction}.fasta", 'fasta')
                for seq_record in msa_record:
                    fname, contig_id = header_parser.match(seq_record.id).groups()
                    cluster_idx = int(cluster_idx)
                    flanking_seq = inter_cluster_sequences[f"{fname}:<-file:contig->:{contig_id}"]
                    flanks.append(flanking_seq[:fname_contigid_to_cluster_dir_to_length[(fname, contig_id)][(cluster_idx, direction)]])
                lengths = Counter(len(s) for s in flanks)
                maxlen = max(lengths)
                if maxlen != 0 and expand_by != 0:
                    bins = np.histogram([len(s) for s in flanks], bins=20)[0]
                    if max(bins)/len(flanks) < 0.99 or lengths.most_common(1)[0][0] == 0:
                        expand_by = lengths.most_common(1)[0][0]
                cluster_space_left = max(cldir_to_len[(cluster_idx, direction)] for cldir_to_len in fname_contigid_to_cluster_dir_to_length.values())
                expand_by = min(
                    expand_by,
                    cluster_space_left
                )
                # expand_by = min(expand_by, 7)
                # expand_by = 10
                # for seq_record in msa_record:
                #     fname, contig_id, cluster_idx = header_parser.match(seq_record.id).groups()
                #     cluster_idx = int(cluster_idx)
                #     expand_by = min(expand_by, fname_contigid_to_cluster_dir_to_length[(fname, contig_id)][(cluster_idx, direction)])
                if expand_by <= 0:
                    continue

                seqs_to_align = {}
                # pprint({key: value.seq[:10] for key, value in inter_cluster_sequences.items()})
                for seq_record in msa_record:
                    fname, contig_id = header_parser.match(seq_record.id).groups()
                    cluster_idx = int(cluster_idx)
                    # print(contig_id)
                    flanking_seq = inter_cluster_sequences[f"{fname}:<-file:contig->:{contig_id}"]
                    # print(flanking_seq.seq)
                    flanking_seq = flanking_seq[:fname_contigid_to_cluster_dir_to_length[(fname, contig_id)][(cluster_idx, direction)]]
                    if len(flanking_seq) >= expand_by:
                        flanking_seq.seq = flanking_seq.seq[:expand_by]
                    seqs_to_align[contig_id] = str(flanking_seq.seq)
                minlen = min(len(s) for s in seqs_to_align.values())
                for k in seqs_to_align:
                    seqs_to_align[k] = seqs_to_align[k][:minlen+indel_cutoff]
                seqlist = list(seqs_to_align.values())
                empty_seqs = [i for i in range(len(seqlist)) if seqlist[i] == ""]
                nonempty_seqs = [s for s in seqlist if s != ""]
                msa_result_temp = []
                if minlen < 50:
                    aligner = pa.msa_aligner()
                    # time.sleep(.25)
                    res = aligner.msa(nonempty_seqs, out_cons=False, out_msa=True)
                    msa_result_temp = res.msa_seq
                else:
                    nonempty_seq_file = f"{cluster_directory}/cluster{cluster_idx}_{direction}_nonempty.fa"
                    SeqIO.write(
                        (SeqIO.SeqRecord(Seq(sequence), id=str(seq_idx)) for seq_idx, sequence in enumerate(nonempty_seqs)),
                        nonempty_seq_file,
                        "fasta")
                    subprocess.check_call(
                        f"muscle -super5 {nonempty_seq_file} -output {nonempty_seq_file}_aligned.fa -threads {cpu_count}", 
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.STDOUT,
                        shell=True)
                    msa_result_temp = [str(record.seq) for record in SeqIO.parse(f"{nonempty_seq_file}_aligned.fa", "fasta")]
                msa_result = []
                idx = 0
                for aligned_seq in msa_result_temp:
                    while idx in empty_seqs:
                        msa_result.append("-"*len(msa_result_temp[0]))
                        idx += 1
                    msa_result.append(aligned_seq)
                    idx += 1
                while idx in empty_seqs:
                    msa_result.append("-"*len(msa_result_temp[0]))
                    idx += 1
                (_, _), ani_cutoff = get_ani_cutoff(msa_result, cutoff)
                msa_seqs = [seq[:ani_cutoff] for seq in msa_result]
                # for i, seq in enumerate(msa_seqs):
                    # print(f">Seq{i}")
                    # print(seq)
                expansion_nucs = max(len(s.replace("-", "")) for s in msa_seqs)
                logger.debug(f"Expanding cluster {cluster_idx} to the {direction} by {expansion_nucs}/{cluster_space_left}")

                # continue
                flanking_seqs = {}
                for idx, seq in enumerate(msa_seqs):
                    contig_id = list(seqs_to_align.keys())[idx]
                    flanking_seqs[contig_id] = Seq(seq)
                for seq_record in msa_record:
                    fname, contig_id = header_parser.match(seq_record.id).groups()
                    cluster_idx = int(cluster_idx)
                    flanking_seq = flanking_seqs[contig_id]
                    if direction == "left":
                        flanking_seq = flanking_seq.reverse_complement()
                        # flanking_seq = "-"*(expand_by - len(flanking_seq)) + flanking_seq
                    # else:
                    #     flanking_seq = flanking_seq + "-"*(expand_by - len(flanking_seq))
                    seq_record.annotations["size"] += len(str(flanking_seq).replace("-", ''))
                    if direction == "right":
                        seq_record.seq = seq_record.seq + flanking_seq
                    else:
                        seq_record.seq = flanking_seq + seq_record.seq
                        seq_record.annotations["start"] -= (len(flanking_seq) - flanking_seq.count("-"))
                        # CHECK IF START > LEN OF SEQ OR SOMETHING
                    adj_cluster_idx, adj_cluster_dir = fname_contigid_to_cluster_dir_to_adjacent_cluster[(fname, contig_id)][(cluster_idx, direction)]
                    fname_contigid_to_cluster_dir_to_length[(fname, contig_id)][(adj_cluster_idx, adj_cluster_dir)] -= len(str(flanking_seq).replace("-", ''))



            # for seq_record in msa_record:
            #     fname, contig_id = header_parser.match(seq_record.id).groups()
                # seq_record.id = contig_id
            msa_record.annotations["cluster"] = cluster_idx
            (record_matches, record_maxmatches), _ = get_match_maxmatch([record.seq for record in msa_record])
            new_matches, new_max_matches = record_matches + new_matches, record_maxmatches + new_max_matches
            new_nucs_aligned += sum(len(record.seq) for record in msa_record)
            write_xmfa_cluster(extended_maf_file, [msa_record], fname_header_to_gcontigidx)
            # maf_writer.write_alignment(msa_record)
    total_length = sum(l for l in fname_contigid_to_length.values())
    logger.debug(f"COVERAGE {(old_nucs_aligned / total_length) * 100:.3f}% --> {(new_nucs_aligned / total_length)*100:.3f}%")
    logger.debug(f"ANI {100*(old_matches/old_max_matches):.2f}% --> {100*(new_matches/new_max_matches):.2f}%")

#%%
def check_maf(maf_file, fname_to_seqrecords):
    contig_to_coords_test = defaultdict(list)
    header_parser = re.compile(r"(.+):<-file:contig->:(.+)")
    with open(maf_file, 'r') as maf_in:
        maf_iterator = MafIterator(maf_in)
        for msa_record in maf_iterator:
            for seq_record in msa_record:
                strand = seq_record.annotations["strand"]
                coord_start = seq_record.annotations["start"]  if seq_record.annotations["strand"] == 1 else  seq_record.annotations["srcSize"] - seq_record.annotations["start"]
                start = seq_record.annotations["start"]
                end = start + seq_record.annotations["size"]
                if strand == -1:
                    coord_end = coord_start - seq_record.annotations["size"]
                    coord_start, coord_end = coord_end, coord_start
                else:
                    coord_end = end
                fname, contig_id = header_parser.match(seq_record.id).groups()
                contig_to_coords_test[(fname, contig_id)].append((coord_start, coord_end, seq_record.annotations["strand"]))
                # if len(seq_record) > 20:
                #     print(seq_record.id.split(":")[0], seq[:15],"...", seq[-15:])
                # else:
                #     print(seq_record.id.split(":")[0], seq)
                block_seq = str(seq_record.seq).lower().replace("-", '').replace("n", ".")
                original_seq = fname_to_seqrecords[fname][contig_id].seq
                if strand == 1:
                    match_pos = re.search(block_seq, str(original_seq), re.IGNORECASE)
                else:
                    original_seq = str(original_seq.reverse_complement())
                    match_pos = re.search(block_seq, original_seq, re.IGNORECASE)
                match_pos = match_pos.start() if match_pos else -1
                if start != match_pos:
                    print(start, "!=", match_pos)
                    print(seq_record.id)
                    print(seq_record.annotations)
                    print(block_seq)
                    return
    for key in contig_to_coords_test:
        contig_to_coords_test[key] = sorted(contig_to_coords_test[key])
    return contig_to_coords_test

def validate_coords(contig_to_coords):
    for contig, coords in contig_to_coords.items():
        prev_end = -1
        for start, end, strand in coords:
            if end <= prev_end:
                print(f"{contig} has invalid coordinates!")
                print(coords)
                break
            prev_end = end

# coords_extended = check_maf("parsnp-extended.maf")
# validate_coords(coords_extended)
# coords_original = check_maf("parsnp.maf")
#%%
# match_score = 4
# mismatch_penalty = -8
# gap_penalty = -10
# validate = True
# genome_dir = "sars"
# sequence_files = glob(f"{genome_dir}/*.fa")
# xmfa_file = f"{genome_dir}/parsnp.xmfa"
# cluster_directory = "clusters/"
# # ref = f"{genome_dir}/ST-11_RT-078_M120.fa"
# ref = f"{genome_dir}/GISAID_9381.fa.ref"
# finalfiles = list(set(sequence_files) - set([ref]))
# temp_directory = "temp/"
# original_maf_file = f"{genome_dir}/parsnp.maf"
# extended_xmfa_file = f"{genome_dir}/parsnp-extended.xmfa"
# fname_contigid_to_length, fname_contigidx_to_header, fname_to_seqrecord = get_sequence_data(
#         ref,
#         finalfiles,
#         index_files=True)
# fname_to_contigid_to_coords, fname_header_to_gcontigidx = xmfa_to_maf(
#         xmfa_file,
#         original_maf_file,
#         fname_contigidx_to_header,
#         fname_contigid_to_length)
# packed_write_result = write_intercluster_regions(finalfiles + [ref], temp_directory, fname_to_contigid_to_coords)
# fname_contigid_to_cluster_dir_to_length, fname_contigid_to_cluster_dir_to_adjacent_cluster = packed_write_result
# cluster_files = glob(f"{temp_directory}/*.fasta")
# clusterdir_expand, clusterdir_len = get_new_extensions(
#         cluster_files,
#         match_score,
#         mismatch_penalty,
#         gap_penalty)
# write_extended_xmfa(
#         original_maf_file,
#         extended_xmfa_file,
#         temp_directory,
#         clusterdir_expand,
#         clusterdir_len,
#         fname_contigid_to_cluster_dir_to_length,
#         fname_contigid_to_cluster_dir_to_adjacent_cluster,
#         fname_header_to_gcontigidx,
#         fname_contigid_to_length)
# parsnp_output = extended_xmfa_file









# fname_contigid_to_length, fname_contigidx_to_header, fname_to_seqrecord = get_sequence_data(sequence_files, index_files=validate)
# fname_to_contigid_to_coords = xmfa_to_maf(xmfa_file, original_maf_file, fname_contigidx_to_header, fname_contigid_to_length)
# packed_write_result = write_intercluster_regions(cluster_directory, fname_to_contigid_to_coords)
# fname_contigid_to_cluster_dir_to_length, fname_contigid_to_cluster_dir_to_adjacent_cluster = packed_write_result
# cluster_files = glob(f"{cluster_directory}/*.fasta")
# clusterdir_expand, clusterdir_len = get_new_extensions(cluster_files)
# write_extended_maf(
#         original_maf_file,
#         extended_maf_file,
#         cluster_directory,
#         clusterdir_expand,
#         clusterdir_len,
#         fname_contigid_to_cluster_dir_to_length,
#         fname_contigid_to_cluster_dir_to_adjacent_cluster)

# if validate:
    # # fname_to_contigid_to_coords = xmfa_to_maf(
    # #         extended_xmfa_file,
    # #         extended_xmfa_file + ".maf",
    # #         fname_contigidx_to_header,
    # #         fname_contigid_to_length)
    # coords_extended = check_maf(extended_xmfa_file + ".maf", fname_to_seqrecord)
    # validate_coords(coords_extended)
    # # coords_original = check_maf(original_maf_file, fname_to_seqrecord)
    # # validate_coords(coords_original)
