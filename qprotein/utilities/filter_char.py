"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/10/11

# Description: 
# ------------------------------------------------------------------------------
"""
from Bio import SeqIO

	
	
def read_save(input_fasta):
	with open(r"C:\Users\bj600\Desktop\db_ncbi.txt", "r") as f:
		db_list = [i.strip() for i in f.readlines()]
	fasta_sequences = SeqIO.parse(open(input_fasta), 'fasta')
	target_ncbi = []
	
	for seq_record in fasta_sequences:
		ncbi_id = seq_record.description.split("_")[1]
		if ncbi_id in db_list:
			target_ncbi.append(ncbi_id)
	for i in target_ncbi:
		print(i)

input_fasta = r"C:\Users\bj600\Desktop\gh11.fasta"
read_save(input_fasta)
	
	