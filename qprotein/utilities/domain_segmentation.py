# -*- coding: utf-8 -*-
# @Time    : 2023/2/24 11:02
# @Author  : Zhixin Dou
"""
Domain segmentation via cath database blastp
"""

import os
from Bio import SeqIO


def read_file(file):
    with open(file, 'r') as f:
        content = f.readlines()
    return content


def read_blastp_output(blastp_output):
    content = read_file(blastp_output)
    blastp_output = []
    for line in content:
        one_output = {}
        line_split = line.rstrip().split('\t')
        one_output['query_name'] = line_split[0]
        one_output['subject_name'] = line_split[1]
        one_output['match_length'] = line_split[3]
        one_output['start_query'] = line_split[6]
        one_output['end_query'] = line_split[7]
        blastp_output.append(one_output)
    print('Step 1: blastp_output have loaded.')
    return blastp_output


def read_fasta(fasta_file):
    fasta = []
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):
        _fasta = [str(seq_record.id), str(seq_record.seq)]
        fasta.append(_fasta)
    print('Step 2: fasta_file have loaded.')
    return fasta


def domain_segmentation(fasta, blastp_output):
    segment_fasta = []
    for record in fasta:
        name = record[0]
        for hit in blastp_output:
            if hit['query_name'] == name:
                seq = record[1][int(hit['start_query']): int(hit['end_query']) + 1]
                _fasta = [name, seq]
                segment_fasta.append(_fasta)
    print('Step 3: Protein domains have been segmented.')
    return segment_fasta


def write_segment_domain(segment_fasta, segment_domain_output_file):
    with open(segment_domain_output_file, 'w') as f:
        for seq_record in segment_fasta:
            f.write('>' + seq_record[0] + '\n')
            f.write(seq_record[1] + '\n')
    print('Step 4: Write output to fasta file')


def parse_cath_sql():
    # TODO: CATH annotation might be useful.
    pass


def main():
    fasta_file = r'D:\subject\active\PyMulstruct\data\GH10\GH10.fasta'
    cath_output = r'D:\subject\active\PyMulstruct\data\GH10\GH10_diamond_filter.out'
    segment_domain_output_file = r'D:\subject\active\PyMulstruct\data\GH10\GH10_segment_domain.fasta'
    blastp_output = read_blastp_output(cath_output)
    fasta = read_fasta(fasta_file)
    segment_fasta = domain_segmentation(fasta, blastp_output)
    write_segment_domain(segment_fasta, segment_domain_output_file)
    print('Domain segmentation finished')


if __name__ == '__main__':
    main()
