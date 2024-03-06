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


def cross_ref(all_temperature, fasta, positive_fasta, negative_fasta):
    filtered_records = []
    high_records = []
    low_records = []
    fasta_sequences = list(SeqIO.parse(open(fasta), 'fasta'))
    for seq_record in fasta_sequences:
        description = seq_record.description
        organism = "_".join(description.split("[")[-1].rstrip("]").split(" ")[:2]).lower()
        for temper in all_temperature:
            if organism == temper[0]:
                seq_record.description = seq_record.description + "_" + temper[2]
                if len(str(seq_record.seq)) > 280:
                    if int(temper[2]) > 50:
                        high_records.append(seq_record)
                    if int(temper[2]) < 40:
                        low_records.append(seq_record)
    print(high_records)
    print(len(low_records))
    # SeqIO.write(high_records, positive_fasta, "fasta")
    # SeqIO.write(low_records, negative_fasta, "fasta")

csv_file = r"D:\subject\active\1-qProtein\data\temperature_db\temperature.csv"
all_temperature = read_temper(csv_file)
fasta_file = r"D:\subject\active\1-qProtein\data\enzymes\GH8\1_preprocessing\5_GH8_rmdup_one_hit_sequence.fasta"
positive_fasta = r"D:\subject\active\1-qProtein\data\enzymes\GH8\1_preprocessing\6_GH48_positive.fasta"
negative_fasta = r"D:\subject\active\1-qProtein\data\enzymes\GH8\1_preprocessing\6_GH48_negative.fasta"
cross_ref(all_temperature, fasta_file, positive_fasta, negative_fasta)
