"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/07/12

# Description: 
# ------------------------------------------------------------------------------
"""
import os
import io

import numpy as np
from tqdm import tqdm

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBIO import PDBIO
import freesasa


def cif2pdb(struct_path):
	struct_name = os.path.splitext(os.path.basename(struct_path))[0]
	parser = MMCIFParser()
	structure = parser.get_structure(struct_name, struct_path)
	pdb_output = io.StringIO()
	pdb_io = PDBIO()
	pdb_io.set_structure(structure)
	pdb_io.save(pdb_output)
	pdb_file_text = io.StringIO(pdb_output.getvalue())
	return pdb_file_text

def get_sasa(structure_path):
	if structure_path.split(".")[-1] == "cif":
		pdb_file_text = cif2pdb(structure_path)
		parser = PDBParser()
		structure = parser.get_structure("pdb", pdb_file_text)
		result, sasa_classes = freesasa.calcBioPDB(structure)
		return result.totalArea(), sasa_classes
	if structure_path.split(".")[-1] == "pdb":
		structure = freesasa.Structure(structure_path)
		result = freesasa.calc(structure)
		sasa_classes = freesasa.classifyResults(result, structure)
		return result.totalArea(), sasa_classes
	
def get_DisulfideBond():
    pass

if __name__ == '__main__':
	structure_dir = r"D:\subject\active\1-qProtein\data\enzymes\cazy_endo-1_4-beta-xylanase\classification\thermophilic"
	file_list = os.listdir(structure_dir)
	polar = []
	for file in tqdm(file_list, desc="Processing"):
		structure_path = os.path.join(structure_dir, file)
		total_area, classes_area = get_sasa(structure_path)
		polar.append(classes_area["Polar"])
	print(np.average(polar))