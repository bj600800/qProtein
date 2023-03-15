#!/usr/bin/env python           
# -*- coding:utf-8 -*-          
# @Filename:    pdb_parser.py      
# @Author:      Eric Dou        
# @Time:        2021/5/27 9:02 

"""Function: Information extraction from .pdb file, i.e. sequence, dss, etc."""
import os

from Bio.PDB import MMCIFParser
from Bio.PDB.DSSP import DSSP
import pdb2pqr

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import PDBParser

def mmcif_dssp(struct_path_list, dssp_path):
    parser = MMCIFParser()

    for struct in struct_path_list:

        mmcif_dict = MMCIF2Dict(struct)
        print(mmcif_dict.keys())
        parser = PDBParser(QUIET=True)
        struct_file_name_pref = os.path.split(struct)[1].split('.')[0]
        struct_name = '_'.join(struct_file_name_pref.split('_')[:6])
        structure = parser.get_structure(struct_name, struct)
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save('example.pdb')
        input()

        struct_file_name_pref = os.path.split(struct)[1].split('.')[0]
        struct_name = '_'.join(struct_file_name_pref.split('_')[:6])
        struct_path = struct

        structure = parser.get_structure(struct_name, struct_path)
        model = structure[0]
        dssp = DSSP(model, struct_path, dssp_path)
        a_key = list(dssp.keys())
        ss = []
        for i in a_key:
            print(dssp[i])
            ss.append(dssp[i][2])




def main():
    struct_dir = r'D:\subject\active\1-qProtein\data\tibet\ident90'
    dssp = r'D:\subject\active\1-qProtein\code\qprotein\structure\tools\dssp\dssp-3.0.0-win32.exe'
    struct_path_list = [os.path.join(struct_dir, i) for i in os.listdir(struct_dir)]
    mmcif_dssp(struct_path_list, dssp)


# main()

from Bio.PDB import PDBParser, MMCIFIO
p = MMCIFParser()
struc = p.get_structure("1446", r"C:\Users\bj600\Desktop\A10_k97_11432_gene_1_1_trembl_A0A1M6UVC4.cif")
io = PDBIO()
io.set_structure(struc)
io.save("file.pdb")