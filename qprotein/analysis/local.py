"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/11/01

# Description: analyze hydrophobic distribution in different local regions
# ------------------------------------------------------------------------------
"""
from pathlib import Path
import csv
import biotite.structure.io as strucio
from biotite.structure import CellList
import numpy as np
from qprotein.utilities import align_map, logger


logger = logger.setup_log(name=__name__)


def get_region(structure, active_archi_res, active_edge_dist, intermediate_edge_dist):
	shell_coord = []
	for shell_res in active_archi_res:
		first_shell_coord = structure[
			(structure.res_id == shell_res) &
			(structure.atom_name == "CA")].coord[0]
		shell_coord.append(first_shell_coord)
	center_coord = np.mean(shell_coord, axis=0)
	cell_list = CellList(atom_array=structure, cell_size=3)  # size relates to the computational efficiency
	active_region = [int(i) for i in list(structure[cell_list.get_atoms(center_coord, radius=active_edge_dist)].res_id)]
	inter_region = [int(i) for i in list(set(structure[cell_list.get_atoms(center_coord, radius=intermediate_edge_dist)].res_id) - set(active_region))]
	surface_region = [int(i) for i in list(set(structure.res_id)-set(active_region)-set(inter_region))]
	
	return active_region, surface_region

def extract_all_cluster_res(cluster_dict):
	all_res = []
	for cluster_info in cluster_dict.values():
		_, residues = cluster_info  # ['326.85', [56,66]]
		all_res.extend(residues)
	return all_res

def analyze_hydrophobic(cluster, active_region, surface_region):
	active_region_hydrophobic = []
	surface_region_hydrophobic = []

	for cluster_id, res_info in cluster.items():
		intersection_first_region_res = list(set(res_info[1]) & set(active_region))
		_first_area = []
		for i in intersection_first_region_res:
			res_list = res_info[1]
			res_id = res_list.index(i)
			res_area = res_info[2][res_id]
			_first_area.append(res_area)
		first_area_sum = sum(_first_area)
		active_region_hydrophobic.append(first_area_sum)

		intersection_second_region_res = list(set(res_info[1]) & set(surface_region))
		_second_area = []
		for i in intersection_second_region_res:
			res_list = res_info[1]
			res_id = res_list.index(i)
			res_area = res_info[2][res_id]
			_second_area.append(res_area)
		second_area_sum = sum(_second_area)
		surface_region_hydrophobic.append(second_area_sum)

	normalized_active_region_hydrophobic = sum(active_region_hydrophobic) / len(active_region)
	normalized_surface_region_hydrophobic = sum(surface_region_hydrophobic) / len(surface_region)
	return normalized_active_region_hydrophobic, normalized_surface_region_hydrophobic


def analyze_aa(structure, active_region, surface_region):
	charge_aa = ["ARG", "LYS", "ASP", "GLU"]
	polar_aa = ["ALA", "GLY", "SER", "THR", "ASN", "GLN", "HIS"]
	hydrophobic_aa = ["CYS", "ILE", "LEU", "VAL", "MET", "PHE", "TRP", "TYR", "PRO"]
	active_aa = []
	surface_aa = []
	overall_aa = list({atom.res_id: atom.res_name for atom in structure}.values())  # using biotite get seq
	
	for res in active_region:
		res_name = structure[structure.res_id == res].res_name[0]
		active_aa.append(res_name)
	for res in surface_region:
		res_name = structure[structure.res_id == res].res_name[0]
		surface_aa.append(res_name)
	
	active_fraction = {}
	distant_fraction = {}
	overall_fraction = {}
	for index, aa_list in enumerate([active_aa, surface_aa, overall_aa]):
		charge_count = 0
		polar_count = 0
		hydrophobic_count = 0
		for aa in aa_list:
			if aa in charge_aa:
				charge_count += 1
			elif aa in polar_aa:
				polar_count += 1
			elif aa in hydrophobic_aa:
				hydrophobic_count += 1
		if index == 0:
			total_count = len(active_aa)
			charge_fraction = charge_count / total_count
			polar_fraction = polar_count / total_count
			hydrophobic_fraction = hydrophobic_count / total_count
			active_fraction["charge"] = charge_fraction
			active_fraction["polar"] = polar_fraction
			active_fraction["hydrophobic"] = hydrophobic_fraction
		if index == 1:
			total_count = len(surface_aa)
			charge_fraction = charge_count / total_count
			polar_fraction = polar_count / total_count
			hydrophobic_fraction = hydrophobic_count / total_count
			distant_fraction["charge"] = charge_fraction
			distant_fraction["polar"] = polar_fraction
			distant_fraction["hydrophobic"] = hydrophobic_fraction
		if index == 2:
			total_count = len(overall_aa)
			charge_fraction = charge_count / total_count
			polar_fraction = polar_count / total_count
			hydrophobic_fraction = hydrophobic_count / total_count
			overall_fraction["charge"] = charge_fraction
			overall_fraction["polar"] = polar_fraction
			overall_fraction["hydrophobic"] = hydrophobic_fraction
	return active_fraction, distant_fraction, overall_fraction
	
	


def save_hydrophobic_feature(hydrophobic_csv_path, data):
	with open(hydrophobic_csv_path, 'w', newline='') as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow(['Protein_name', 'Active', 'Distant'])
		for name, values in data.items():
			writer.writerow([name, values[0], values[1]])
	logger.info(f"Local hydrophobic feature saved: {hydrophobic_csv_path}")

def save_aa_feature(aa_csv_path, data):
	with open(aa_csv_path, 'w', newline='') as csvfile:
		writer = csv.writer(csvfile)
		writer.writerow(['Protein_name', 'active_charge', 'active_polar', 'active_hydrophobic',
		                 'distant_charge', 'distant_polar', 'distant_hydrophobic',
		                 'overall_charge', 'overall_polar', 'overall_hydrophobic'])
		for name, region in data.items():
			writer.writerow([name, str(region["active"]["charge"]), str(region["active"]["polar"]), str(region["active"]["hydrophobic"]),
			                 str(region["distant"]["charge"]), str(region["distant"]["polar"]), str(region["distant"]["hydrophobic"]),
			                 str(region["overall"]["charge"]), str(region["overall"]["polar"]), str(region["overall"]["hydrophobic"])])
	logger.info(f"Local amino acid feature saved: {aa_csv_path}")

def run(template_name, alignment_fasta, template_active_architecture, clusters, pdb_list,
        hydrophobic_csv_path, aa_csv_path, active_edge_dist, intermediate_edge_dist):
	map_dict = align_map.run(alignment_fasta)
	aligned_active_architecture = [map_dict[template_name][int(i)] for i in template_active_architecture]
	hydrophobic_region_dict = {}
	aa_region_dict = {}
	for structure_path in pdb_list:
		structure_name = Path(structure_path).stem
		structure = strucio.load_structure(structure_path)
		res_map = map_dict.get(structure_name)
		if res_map:
			reversed_res_map = {v:k for k, v in res_map.items()}
			active_archi_res = [reversed_res_map.get(i) for i in aligned_active_architecture]
			active_archi_res = [i for i in active_archi_res if i is not None]  # clean None
			active_region, surface_region = get_region(structure, active_archi_res, active_edge_dist, intermediate_edge_dist)
			hydrophobic_region_dict[structure_name] = analyze_hydrophobic(clusters[structure_name], active_region, surface_region)
			aa_ret = analyze_aa(structure, active_region, surface_region)
			aa_dict = {}
			aa_dict["active"] = aa_ret[0]
			aa_dict["distant"] = aa_ret[1]
			aa_dict["overall"] = aa_ret[2]
			aa_region_dict[structure_name] = aa_dict
	save_hydrophobic_feature(hydrophobic_csv_path, hydrophobic_region_dict)
	save_aa_feature(aa_csv_path, aa_region_dict)

if __name__ == '__main__':
	template_name = r"A5H0S3"
