from functools import reduce
import operator
import os
import re
from tqdm import tqdm
import argparse
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import csv


# python verify.py parsnp.xmfa ./ ./test/ -m -f

# I made them positional arguments so that there's no need of dashes
# also I feel like we could just use the file names in xmfa? that way we
# only need to provide path to all the fna files; not sure

parser = argparse.ArgumentParser()
parser.add_argument('xmfa', type=str, help="the path to the xmfa file you are trying to verify")
parser.add_argument('ref', type=str, help="the path to the referrence fna file")
parser.add_argument('fna', type=str, help="the path to all the fna files")
parser.add_argument('-m', '--maf', action='store_true', help="exports a .maf file translated from the xmfa file")
parser.add_argument('-f', '--find_actual', action='store_true', 
    help="find the actual coordinates of the misplaced alignments. Slowing it down significantly")
args = parser.parse_args()

xmfa_path = args.xmfa
ref_path = args.ref
fna_path = args.fna
maf_flag = args.maf
find_flag = args.find_actual

def compare_with_dashes(str1, str2):
    # ignores the dashes when comparing
    if str1 == str2:
        return True
    if len(str1) != len(str2):
        return False
    else:
        return all(c1 == c2 for (c1, c2) in filter(lambda pair: '-' not in pair, zip(str1, str2))) or \
        all(c1 == c2 for (c1, c2) in zip(filter(lambda x: x!='-', str1),(filter(lambda x: x!='-', str2))))
        # deals with gaps


def reverse_complement(str1):
    return str(Seq(str1).reverse_complement())


with open(xmfa_path) as xmfa:
    line = xmfa.readline()
    line = xmfa.readline()
    seqNum = int(line.split()[1])
    seqs = {}
    for i in tqdm(range(seqNum)):
        line = xmfa.readline()
        line = xmfa.readline()
        seqName = line.split()[1]
        if seqName[-4:] == ".ref":
            seqName = seqName[:-4]
        seqs[i+1] = seqName
        line = xmfa.readline()
        line = xmfa.readline()

    line = xmfa.readline()
    intervalCount = int(line.split()[1])

    # Header parsing over

    seqVerify = [] # a dictionary 
    for seq in range(seqNum+1):
        seqVerify.append ([])

    line = xmfa.readline()

    with tqdm(total=intervalCount) as pbar:
        while line:
            if line.find(" + ") != -1:
                alignment = re.split("-|:p| cluster| s|:|\s", line[1:])
            else:
                alignment = re.split("(?<!\s)-|:p| cluster| s|:|\s", line[1:])
                # here regex doesn't match dashes that are preceded by spaces
            # Here the alignments are in order:
            # [seqeunce number, starting coord, end coord, ...
            # + or -, cluster number, contig number, coord in contig]
            # TODO parse the negative coords correctly
            line = xmfa.readline()
                # here only forward alignments are used
                # here we store the contig and coord relative to contig.
            seqVerify[int(alignment[0])].append((int(alignment[2])-int(alignment[1]),alignment[3],int(alignment[5]), int(alignment[6]), line[:20]))

            # notice that here only the first 20 are taken
            # get to next alignment header
            while line and (initial := line[0]) != '>':
                if initial == '=':
                    pbar.update(1)
                line = xmfa.readline()

# parsing xmfa done.

now = datetime.now()

current_time = now.strftime("%Y-%m-%d-%H%M%S")

contig_size = {}

dnas = []
dnas.append([])
for seq in range(1,seqNum+1):
    if seq == 1:
        path = ref_path + seqs[seq]
    else:
        path = fna_path + seqs[seq]
    dnas.append ([record.seq for record in SeqIO.parse(path, "fasta")]) #input file 

counter = 0
csvdata=[]
with open(current_time+".txt", "x") as f:
    for idx in range(len(seqVerify[1])): #idx is the index of alignment
        dicts=[]
        dicts.append({})
        for seq in range(1,seqNum+1):
            dna=dnas[seq]
            contig_size[seq] = {i+1: len(contig) for i, contig in enumerate(dna)} 
            # here we store the size of each contig (index starting with 1)
            # to prepare for the .maf file
            coords = seqVerify[seq][idx]
            alignment_length = coords[0]
            strand = coords[1]
            contig = coords[2]
            target = coords[3]
            xmfa_seq = coords[4]
            #for alignment_length, strand, contig, target, xmfa_seq in coords:
            dict={}
            dict["sequence"] = seq
            dict["file name"] = seqs[seq]
            dict["strand"] = strand
            dict["contig number"] = contig
            dict["contig position"] = target 
            seq_length = contig_size[seq][contig]
            length = len(xmfa_seq.strip()) 
            if contig == 1 and target ==1:
                target = 0 
            if strand == '+':
                fna_seq=dna[contig-1][target:target+length].lower()
            if strand == '-':
                fna_seq=reverse_complement(dna[contig-1].lower()[target-length:target])
            if len(fna_seq) == len(xmfa_seq) and not compare_with_dashes(fna_seq, xmfa_seq.lower()):
                f.write("sequence: " + str(seq) + "\n")
                f.write("file name: " + seqs[seq] + "\n")
                f.write("strand: " + str(strand) + "\n")
                f.write("position in xmfa: s" + str(contig) + ":p" + str(target) + "\n")
                dict["difference between actual position and xmfa output"] = None
                actual_pos = None
                if find_flag:
                    if strand == '+':
                        for cont, s in enumerate(dna):
                            if (pos:=s.lower().find(xmfa_seq.lower())) != (-1):
                                actual_pos = (cont+1, pos+1)
                    else:
                        if (pos:=(dna[contig-1].lower().find(reverse_complement(xmfa_seq.lower())))) != (-1):
                            actual_pos = (contig, pos)
                            # find in reverse conplement for reverse strands 
                    if actual_pos:
                        f.write("actual position: s" + str(actual_pos[0]) + ":p" + str(actual_pos[1]) + '\n')
                        dict["actual contig number"] = actual_pos[0]
                        dict["actual contig position"] = actual_pos[1]
                        f.write(f"Seq length:{seq_length}, target:{target}, alignment length:{alignment_length}, xmfa length:{length} \n")
                        if strand == '-':
                            f.write(f"{(target-length)-actual_pos[1]}\n")
                            if actual_pos[0] == contig:
                                dict["difference between actual position and xmfa output"] = (target-length)-actual_pos[1]
                        else:
                            f.write(f"{target-1-actual_pos[1]}\n")
                            if actual_pos[0] == contig:
                                dict["difference between actual position and xmfa output"] = target-1-actual_pos[1]
                f.write("fna: " + str(fna_seq) + "\n")
                dict["fna"] = str(fna_seq)
                f.write("xmfa: " + xmfa_seq.lower() + "\n")
                dict["xmfa"] = xmfa_seq.lower()
                f.write("----" + "\n")
            else:
                dict["actual contig number"] = None
                dict["actual contig position"] = None
                dict["difference between actual position and xmfa output"] = 0
                dict["fna"] = str(fna_seq)
                dict["xmfa"] = xmfa_seq.lower()
            dicts.append(dict)
        for dict in dicts:
            if dict != {}:
                csvdata.append(dict)

#print(counter)

#print(csvdata)

fieldnames = csvdata[0].keys()
csv_file_path = current_time+".csv"
with open(csv_file_path, 'w', newline='') as file:
    writer = csv.DictWriter(file, fieldnames=fieldnames)

    # Write column headers
    writer.writeheader()

    # Write the data rows
    writer.writerows(csvdata)

if maf_flag:
    data = []
    with open(xmfa_path, "r+") as xmfa:
        for line in xmfa:
            if line[0] == '>' and line[1] != " ":
                data.append("> " + line[1:])
            else:
                data.append(line)

    with open("temp_xmfa", 'x') as tmp:
        for d in data:
            tmp.write(d)

    alignments = AlignIO.parse("temp_xmfa", 'mauve')
    # starting to create new MultipleSequenceAlignment objects
    # here only start, strand and srcSize is needed for maf to format
    # size of alignment is calculated by Biopython
    # Only thing we don't just get from original id is file name and contig size.
    new_msas = []
    for block, aln in enumerate(alignments):
        msa = []
        for index, seq in enumerate(aln):
            og_id = re.split(" s|:p|/", seq.id)
            # [cluster, contig index, pos relative to contig, absolute pos]
            new_seq = SeqRecord(
                seq.seq,
                name = seq.name,
                id = seqs[index+1].split(".")[0] + ":" + og_id[1],
                annotations={
                    "start": og_id[2],
                    "strand": seq.annotations["strand"],
                    "srcSize": contig_size[index+1][int(og_id[1])]
                }
            )
            msa.append(new_seq)
        align = MultipleSeqAlignment(msa)
        new_msas.append(align)
    
    with open(current_time+".maf", "x") as file:
        maf = AlignIO.MafIO.MafWriter(file)
        for align in new_msas:
            maf.write_alignment(align)

    os.remove("temp_xmfa")