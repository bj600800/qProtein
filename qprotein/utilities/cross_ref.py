"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/06/28

# Description: temperature annotation via cross-reference
# ------------------------------------------------------------------------------
"""
import csv
from Bio import SeqIO


def read_temper(file):
    all_temperature = []
    with open(file, "r") as f:
        csv_reader = csv.reader(f)
        for row in csv_reader:
            all_temperature.append(row)
    return all_temperature


def cross_ref(all_temperature, fasta, output_fasta):
    filtered_records = []
    high_records = []
    fasta_sequences = list(SeqIO.parse(open(fasta), 'fasta'))
    for seq_record in fasta_sequences:
        description = seq_record.description
        organism = "_".join(description.split("[")[-1].rstrip("]").split(" ")[:2]).lower()
        for temper in all_temperature:
            if organism == temper[0]:
                seq_record.description = seq_record.description + "_" + temper[2]
                if int(temper[2]) > 45:
                    high_records.append(seq_record.description)
                filtered_records.append(seq_record)
    print(len(high_records))
    SeqIO.write(filtered_records, output_fasta, "fasta")

csv_file = r"D:\subject\active\1-qProtein\data\bacteria_temperature.csv"
all_temperature = read_temper(csv_file)
fasta_file = r"D:\subject\active\1-qProtein\data\enzymes\cazy_endo-1_4-beta-xylanase\3_crawl_sequences_GH11_sequences.fasta"
output_fasta = r"D:\subject\active\1-qProtein\data\enzymes\cazy_endo-1_4-beta-xylanase\3_temperature_sequences_GH11_sequences.fasta"
cross_ref(all_temperature, fasta_file, output_fasta)
