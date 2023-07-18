"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/07/11

# Description:
The implementation of the following detect_disulfide_bonds function is based on the code by Sphinx-Gallery
Original code source: https://www.biotite-python.org/examples/gallery/structure/disulfide_bonds.html
The modified code is used to detect disulfide bonds in protein structures and returns the frequency of disulfide bonds.
# ------------------------------------------------------------------------------
"""
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as pdb
from tqdm import tqdm
import os

def detect_disulfide_bonds(structure, distance=2.05, distance_tol=0.05, dihedral=90, dihedral_tol=15):
	disulfide_bonds = []
	sulfide_mask = (structure.res_name == "CYS") & \
				   (structure.atom_name == "SG")
	cell_list = struc.CellList(
			structure,
			cell_size=distance + distance_tol,
			selection=sulfide_mask
			)
	for sulfide_i in np.where(sulfide_mask)[0]:
		potential_bond_partner_indices = cell_list.get_atoms_in_cells(
				coord=structure.coord[sulfide_i]
				)
		for sulfide_j in potential_bond_partner_indices:
			if sulfide_i == sulfide_j:
				continue
			sg1 = structure[sulfide_i]
			sg2 = structure[sulfide_j]
			cb1 = structure[
				(structure.chain_id == sg1.chain_id) &
				(structure.res_id == sg1.res_id) &
				(structure.atom_name == "CB")
				]
			cb2 = structure[
				(structure.chain_id == sg2.chain_id) &
				(structure.res_id == sg2.res_id) &
				(structure.atom_name == "CB")
				]
			bond_dist = struc.distance(sg1, sg2)
			bond_dihed = np.abs(np.rad2deg(struc.dihedral(cb1, sg1, sg2, cb2)))
			if (distance - distance_tol < bond_dist < distance + distance_tol and
					dihedral - dihedral_tol < bond_dihed < dihedral + dihedral_tol):
				bond_tuple = sorted((sulfide_i, sulfide_j))
				if bond_tuple not in disulfide_bonds:
					disulfide_bonds.append(bond_tuple)
	return np.array(disulfide_bonds, dtype=int)


def read_structure(struct_path):
	structure = ''
	if struct_path.split('.')[-1] == 'cif':
		cif_file = pdbx.PDBxFile.read(struct_path)
		structure = pdbx.get_structure(cif_file)[0]
	if struct_path.split('.')[-1] == 'pdb':
		pdb_file = pdb.PDBFile.read(struct_path)
		structure = pdb_file.get_structure()[0]
	return structure


def run(structure_dir):
	bond_num = []
	structure_list = os.listdir(structure_dir)
	for file in tqdm(structure_list, desc="Progress"):
		structure_file = os.path.join(structure_dir, file)
		structure = read_structure(structure_file)
		disulfide_bonds = detect_disulfide_bonds(structure)
		bond_num.append(len(disulfide_bonds))
	return bond_num

# structure_dir = r"D:\subject\active\1-qProtein\data\enzymes\cazy_endo-1_4-beta-xylanase\classification\thermophilic"
structure_dir = r"D:\subject\active\1-qProtein\data\enzymes\cazy_endo-1_4-beta-xylanase\classification\mesophilic"

bonds_num = run(structure_dir)
frequency = sum(bonds_num)/len(bonds_num)
print("Disulfide bond frequency:", "{:.3f}".format(frequency))
