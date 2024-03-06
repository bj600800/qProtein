"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/10/11

# Description: 
# ------------------------------------------------------------------------------
"""
import os

from Bio import SeqIO

one_hit = []
dele_acc = []

dir_path = r"C:\Users\bj600\Desktop\GH48"
files = os.listdir(dir_path)
for file in files:
	file_path = os.path.join(dir_path, file)
	with open(file_path, "r") as f:
		content = f.readlines()
		for i in content:
			hit_tag_str = '_'.join(i.split("\t")[1:])
			if hit_tag_str not in one_hit:
				one_hit.append(hit_tag_str)
			else:
				dele_acc.append(i.split("\t")[0])

print("delete num: ", len(dele_acc))
out_fasta = []
fasta_sequences = SeqIO.parse(open(r"D:\subject\active\1-qProtein\data\enzymes\GH48\1_preprocessing\3_GH48_sequence.fasta"), 'fasta')
for seq_record in fasta_sequences:
	if seq_record.description.split(" ")[0] not in dele_acc:
		if len(seq_record.seq) > 600:
			out_fasta.append(seq_record)

SeqIO.write(out_fasta, r"D:\subject\active\1-qProtein\data\enzymes\GH48\1_preprocessing\4_GH48_one_hit_sequence.fasta", "fasta")