#!/usr/bin/env python           
# -*- coding:utf-8 -*-          
# @Filename:    pdb_parser.py      
# @Author:      Eric Dou        
# @Time:        2021/5/27 9:02 

"""Function: Information extraction from .pdb file, i.e. sequence, dss, etc."""
import os

from Bio.PDB import MMCIFParser
from Bio.PDB.DSSP import DSSP


def mmcif_dssp(struct_path_list, ss_output_file):
    parser = MMCIFParser()
    with open(ss_output_file, 'w') as wt:
        for struct in struct_path_list:
            struct_file_name_pref = os.path.split(struct)[1].split('.')[0]
            struct_name = '_'.join(struct_file_name_pref.split('_')[:6])
            struct_path = struct
            # if struct_name in acc_id_list or struct_file_name_pref == '5eba':
            print(struct_name)
            structure = parser.get_structure(struct_name, struct_path)
            model = structure[0]
            dssp = DSSP(model, struct_path,
                        dssp=r'D:\subject\active\1-qProtein\code\qprotein\structure\dssp\dssp-3.0.0-win32.exe'
                        )
            a_key = list(dssp.keys())
            ss = []
            for i in a_key:
                print(dssp[i])
                ss.append(dssp[i][2])
            secondary_struct = ''.join(ss).replace('G', 'H').replace('I', 'H') \
                .replace('T', 'c').replace('S', 'c').replace('B', 'c').replace('-', 'c')

            wt.write('>' + struct_name.lstrip() + '\n')
            wt.write(secondary_struct + '\n')


def main():
    struct_dir = r'D:\subject\active\1-qProtein\data\tibet\structure'
    ss_output_file = r'D:\subject\active\1-qProtein\data\tibet\ss_output.txt'
    struct_path_list = [os.path.join(struct_dir, i) for i in os.listdir(struct_dir)]
    mmcif_dssp(struct_path_list, ss_output_file)


main()
