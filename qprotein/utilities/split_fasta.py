"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/04/15

# Description: Split fasta from one file into a directory
# ------------------------------------------------------------------------------
"""
import os
import argparse

from Bio.SeqIO import parse

parser = argparse.ArgumentParser(description='Split fasta sequences within a file to multiple.')

parser.add_argument('--fasta_path', required=True, help='target fasta file to split')
parser.add_argument('--dir_path', required=True, help='output dir for segmented fasta files')

args = parser.parse_args()


def make_dir(dir_path):
    if not os.path.exists(dir_path):
        print(f'Making directory in {dir_path}')
        os.makedirs(os.path.join(dir_path))


def read_fasta(fasta_path):
    fasta = []
    print('Reading fasta file.')
    for i in parse(open(fasta_path), 'fasta'):
        fasta.append((str(i.id), str(i.seq).replace('*', '')))
    return fasta


def split_fasta(fasta, dir_path):
    make_dir(dir_path)
    print('Splitting sequences into fasta files.')
    for name, seq in fasta:
        name = name.replace('|', '_')
        with open(os.path.join(dir_path, name + '.fasta'), 'w') as f:
            f.write('>' + name + '\n')
            f.write(seq + '\n')


def main():
    fasta_path = args.fasta_path
    dir_path = args.dir_path
    fasta = read_fasta(fasta_path)
    split_fasta(fasta, dir_path)


if __name__ == '__main__':
    main()
