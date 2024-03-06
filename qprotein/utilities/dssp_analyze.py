"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/10/25

# Description: 
# ------------------------------------------------------------------------------
"""
import os
from tqdm import tqdm
import numpy as np
import warnings
from Bio import BiopythonWarning
warnings.simplefilter("ignore", BiopythonWarning)

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP


def get_dssp_dat(struct_path):
    struct_file_name = os.path.split(struct_path)[1].split('.')[0]
    parser = PDBParser()
    model = parser.get_structure(struct_file_name, struct_path)[0]
    # dssp = DSSP(model, struct_path, dssp="mkdssp")
    dssp = DSSP(model, struct_path, dssp=r"D:\subject\active\1-qProtein\tools\dssp\dssp-3.0.0-win32.exe")
    dssp_dat = [residue_dat for residue_dat in dssp]
    return dssp_dat


def get_second_struct(struct_path):
    dssp_dat = get_dssp_dat(struct_path)
    second_struct = [residue_dat[2].replace('-', 'C') for residue_dat in dssp_dat]
    ss = ''.join(second_struct).replace('G', 'H').replace('I', 'H').replace('B', 'E').replace('S', 'C').replace('T', 'C')
    return ss


def get_ss_content(ss):
    ss_content = {}
    total_length = len(ss)
    ss_content['helix'] = ss.count('H')/total_length
    ss_content['sheet'] = ss.count('E')/total_length
    ss_content['loop'] = ss.count('C') / total_length
    # ss_content['turn'] = ss.count('T')/total_length
    return ss_content


def count_ss_lengths(ss):
    ss = ss.replace('T', 'C')
    lengths = []
    temp_str = ""
    prev_char = ""

    for char in ss:
        if char == prev_char:
            temp_str += char
        else:
            if prev_char:
                lengths.append((prev_char, len(temp_str)))
            temp_str = char

        prev_char = char

    # 处理最后一个字符串
    lengths.append((prev_char, len(temp_str)))

    return lengths


def analyze(struct_path):
    ss = get_second_struct(struct_path)
    ss_content = get_ss_content(ss)
    length = count_ss_lengths(ss)
    return ss_content



if __name__ == '__main__':
	import csv
	from tqdm import tqdm
	root_path = r"D:\subject\active\1-qProtein\data\enzymes\GH10\2_StrucMapping"
	ret_path = r"D:\subject\active\1-qProtein\data\enzymes\GH10\3_StrucAnalyzing\10_ss_content.csv"
	structure_dir = ("negative", "positive")
	structure_dir = [os.path.join(root_path, i) for i in structure_dir]
	structure_path = [os.path.join(one_dir, file) for one_dir in structure_dir for file in os.listdir(one_dir)]
	ss_content = []
	for path in tqdm(structure_path):
		out = analyze(path)
		ss_content.append(out.values())
	with open(ret_path, "w", newline='') as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow(["helix", "sheet", "loop"])
		for i in ss_content:
			writer.writerow(i)