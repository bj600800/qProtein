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
    sequence_length_threshold = 200
    fasta_sequences = SeqIO.parse(open(input_fasta), 'fasta')
    with open(output_fasta, "w") as wf:
        for seq_record in fasta_sequences:
            sequence = seq_record.seq
            sequence_length = len(sequence)
            if sequence_length >= sequence_length_threshold:
                SeqIO.write(seq_record, wf, "fasta")


def main():
    input_fasta = r"D:\subject\active\1-qProtein\data\enzymes\endo-1_4-beta-xylanase\GH11_sequences.fasta"
    output_fasta = r"D:\subject\active\1-qProtein\data\enzymes\endo-1_4-beta-xylanase\GH11_sequences_filtered.fasta"
    read_save(input_fasta, output_fasta)


if __name__ == '__main__':
    main()