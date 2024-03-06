"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/07/03

# Description: 
# ------------------------------------------------------------------------------
"""
import os.path

from Bio import SeqIO
from Bio.Seq import Seq


def get_domain(input_fasta, output_fasta):
    begain_posi, end_posi = (3099, 7274)
    domain_records = []
    fasta_sequences = list(SeqIO.parse(open(input_fasta), 'fasta'))
    for seq_record in fasta_sequences:
        seq_record.seq = Seq(str(seq_record.seq[begain_posi:end_posi + 1]).replace("-", ""))
        domain_records.append(seq_record)
    SeqIO.write(domain_records, output_fasta, "fasta")


def main():
    input_fasta = r"D:\subject\active\1-qProtein\data\enzymes\GH8\1_preprocessing\8_GH8_align_merge_temperature.fasta"
    output_fasta = r"D:\subject\active\1-qProtein\data\enzymes\GH8\1_preprocessing\9_GH8_cut_align_merge_temperature.fasta"
    get_domain(input_fasta=input_fasta, output_fasta=output_fasta)


if __name__ == '__main__':
    main()