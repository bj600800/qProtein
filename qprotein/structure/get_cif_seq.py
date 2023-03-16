"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/14

# Description: 
# ------------------------------------------------------------------------------
"""
import os

from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def get_seq(cif):
    mmcif_dict = MMCIF2Dict(cif)

    seq = mmcif_dict['_struct_ref.pdbx_seq_one_letter_code'][0].replace('\n', '')
    # if os.path.split(cif)[1] == 'A1_k97_5727_gene_1_1_trembl_A0A6J7I5M5.cif':
    #     print(seq)
    return seq


def write_fasta(cif_dir, output_fasta):
    files = os.listdir(cif_dir)
    with open(output_fasta, 'w') as wf:
        for file in files:
            cif = os.path.join(cif_dir, file)
            seq = get_seq(cif)
            name = os.path.split(file)[1].split('.')[0]
            wf.write('>' + name + '\n')
            wf.write(seq + '\n')


if __name__ == '__main__':
    cif_dir = r'D:\subject\active\1-qProtein\data\manure\ident90'
    output_fasta = r'D:\subject\active\1-qProtein\data\manure\ident90_from_cif.fasta'
    write_fasta(cif_dir, output_fasta)
