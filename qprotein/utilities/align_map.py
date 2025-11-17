"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/09/19

# Description: Create sequence map for local analysis
# ------------------------------------------------------------------------------
"""
from Bio import SeqIO


def get_seq_map(alignment):
    seq_dict = {}
    with open(alignment, 'r') as file:
        lines = [line for line in file if "#Total CPU time is" not in line]
    with open(alignment, 'w') as file:
        file.writelines(lines)

    for record in SeqIO.parse(alignment, "fasta"):
        name = record.description.split(".pdb")[0].lstrip('/')
        seq_dict[name] = (str(record.seq), str(record.seq))

    return seq_dict


def get_id_map(seq_dict):
	id_map = {}
	for name, seq_tuple in seq_dict.items():
		alignment_dict = {}
		original_index = 0
		aligned_index = 0
		original_sequence = seq_tuple[0]
		aligned_sequence = seq_tuple[1]
		while aligned_index < len(aligned_sequence):
			if aligned_sequence[aligned_index] != "-":
				alignment_dict[original_index + 1] = aligned_index + 1
				original_index += 1
		
			aligned_index += 1
		
			if original_index == len(original_sequence):
				break
		id_map[name] = alignment_dict
	return id_map


def run(alignment_fasta):
	seq_dict = get_seq_map(alignment_fasta)
	id_map = get_id_map(seq_dict)
	return id_map


