"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/06/23

# Description: 
# ------------------------------------------------------------------------------
"""
import os

from Bio import SeqIO


def read_save(input_fasta, output_fasta):
    fasta_sequences = SeqIO.parse(open(input_fasta), 'fasta')
    name = set()
    seq_ret = []
    for seq_record in fasta_sequences:
        organism = "_".join(seq_record.description.split("[")[-1].rstrip("]").split(" ")[:2]).lower()
        if organism not in name:
            seq_ret.append(seq_record)
        name.add(organism)
    
    SeqIO.write(seq_ret, output_fasta, "fasta")
    
        
        


def main():
    input_fasta = r"D:\subject\active\1-qProtein\data\enzymes\GH13_5\1_preprocessing\4_temperature_GH13_5_sequences.fasta"
    output_fasta = r"D:\subject\active\1-qProtein\data\enzymes\GH13_5\1_preprocessing\5_negative_GH13_5_sequences.fasta"
    read_save(input_fasta, output_fasta)


if __name__ == '__main__':
    main()