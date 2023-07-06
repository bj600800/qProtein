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

def get_position(fasta):
    fasta_sequences = list(SeqIO.parse(open(fasta), 'fasta'))
    first_position = []
    last_position = []
    for seq_record in fasta_sequences:
        if "temperature" in seq_record.description:
            print(seq_record.description)
            for position, resi in enumerate(seq_record.seq):
                if resi != "-":
                    first = position
                    first_position.append(first)
                    print("First position: ", first)
                    break
            for position, resi in enumerate(reversed(seq_record.seq)):
                if resi != "-":
                    last = len(seq_record.seq) - position
                    last_position.append(last)
                    print("last position: ", last)
                    break
    return min(first_position), max(last_position)


def get_domain(positions, input_fasta, output_fasta):
    begain_posi, end_posi = positions
    domain_records = []
    fasta_sequences = list(SeqIO.parse(open(input_fasta), 'fasta'))
    for seq_record in fasta_sequences:
        seq_record.seq = Seq(str(seq_record.seq[begain_posi:end_posi + 1]).replace("-", ""))
        domain_records.append(seq_record)
    SeqIO.write(domain_records, output_fasta, "fasta")


def main():
    input_fasta = r"D:\subject\active\1-qProtein\data\enzymes\ncbi_endo-1_4-beta-xylanase\5_merged_thermo_ncbi_aligned.fasta"
    output_fasta = os.path.join(os.path.dirname(input_fasta), "6_temperature_domain.fasta")
    positions = get_position(input_fasta)
    get_domain(positions=positions, input_fasta=input_fasta, output_fasta=output_fasta)


if __name__ == '__main__':
    main()