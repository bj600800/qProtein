"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/11/01

# Description: 
# ------------------------------------------------------------------------------
"""
import os

import biotite.structure.io as strucio
from biotite.structure import CellList
import numpy as np
from tqdm import tqdm
from qprotein.utilities import align_map
from qprotein.characterization import hydrophobic
from itertools import zip_longest
import csv

def get_region(structure, active_archi_res):
	# get coordinate
	shell_coord = []
	for shell_res in active_archi_res:
		first_shell_coord = structure[
			(structure.res_id == shell_res) &
			(structure.atom_name == "CA")].coord[0]
		shell_coord.append(first_shell_coord)
	center_coord = np.mean(shell_coord, axis=0)
	cell_list = CellList(atom_array=structure, cell_size=3)
	active_region = list(structure[cell_list.get_atoms(center_coord, radius=12)].res_id)
	inter_region = list(set(structure[cell_list.get_atoms(center_coord, radius=15)].res_id) - set(structure[cell_list.get_atoms(center_coord, radius=10)].res_id))
	surface_region = list(set(structure.res_id)-set(active_region)-set(inter_region))
	
	return active_region, surface_region

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
	surface_fraction = {}
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
			surface_fraction["charge"] = charge_fraction
			surface_fraction["polar"] = polar_fraction
			surface_fraction["hydrophobic"] = hydrophobic_fraction
		if index == 2:
			total_count = len(overall_aa)
			charge_fraction = charge_count / total_count
			polar_fraction = polar_count / total_count
			hydrophobic_fraction = hydrophobic_count / total_count
			overall_fraction["charge"] = charge_fraction
			overall_fraction["polar"] = polar_fraction
			overall_fraction["hydrophobic"] = hydrophobic_fraction
	return active_fraction, surface_fraction, overall_fraction
	
	
def analyze_hydrophobic(structure, active_region, surface_region):
	hydrophobic_out = hydrophobic.analyze(structure)["cluster_res"]
	active_region_hydrophobic = []
	surface_region_hydrophobic = []
	for out in hydrophobic_out:
		intersection_first_region_res = list(set(out.keys()) & set(active_region))
		active_region_hydrophobic.append(sum([float(out[i]) for i in intersection_first_region_res]))
		
		intersection_second_region_res = list(set(out.keys()) & set(surface_region))
		surface_region_hydrophobic.append(sum([float(out[i]) for i in intersection_second_region_res]))
	normalized_active_region_hydrophobic = sum(active_region_hydrophobic) / len(active_region)
	normalized_surface_region_hydrophobic = sum(surface_region_hydrophobic)/len(surface_region)

	return normalized_active_region_hydrophobic, normalized_surface_region_hydrophobic


def write2csv(csv_path, data):
	max_length = max(len(lst) for lst in data)
	transposed_data = list(zip_longest(*data, fillvalue=None))
	
	with open(csv_path, 'w', newline='') as csvfile:
		writer = csv.writer(csvfile)
		# 写入数据
		for i, row in enumerate(transposed_data, start=1):
			writer.writerow([i] + list(row[:max_length]))


def run(template_name, template_shell_file, origin_fasta_file, align_fasta_file, structure_dir):
	with open(template_shell_file, "r") as f:
		template_first_shell = [i.rstrip() for i in f.readlines()]
	map_dict = align_map.run(origin_fasta_file, align_fasta_file)
	align_first_shell = [map_dict[template_name][int(i)] for i in template_first_shell]
	all_active_area = []
	all_surface_area = []
	
	all_active_charge_fraction = []
	all_surface_charge_fraction = []
	all_overall_charge_fraction = []
	
	all_active_polar_fraction = []
	all_surface_polar_fraction = []
	all_overall_polar_fraction = []
	
	all_active_hydrophobic_fraction = []
	all_surface_hydrophobic_fraction = []
	all_overall_hydrophobic_fraction = []
	
	for one_dir in structure_dir:
		active_hydrophobic = []
		surface_hydrophobic = []
		
		active_charge_fraction = []
		active_polar_fraction = []
		active_hydrophobic_fraction = []
		
		surface_charge_fraction = []
		surface_polar_fraction = []
		surface_hydrophobic_fraction = []
		
		overall_charge_fraction = []
		overall_polar_fraction = []
		overall_hydrophobic_fraction = []
		
		for file in os.listdir(one_dir):
			structure_path = os.path.join(one_dir, file)
			# print(structure_path)
			structure_name = os.path.splitext(file)[0]
			structure = strucio.load_structure(structure_path)
			res_map = map_dict.get(structure_name)
			if res_map:
				reversed_res_map = {v:k for k, v in res_map.items()}
				active_archi_res = [reversed_res_map.get(i) for i in align_first_shell]
				active_archi_res = [i for i in active_archi_res if i is not None]  # clean None
				active_region, surface_region = get_region(structure, active_archi_res)
				# print("+".join([str(i) for i in active_region]))
				# print("+".join([str(i) for i in surface_region]))
				# print()
				normalized_active_region_hydrophobic, normalized_surface_region_hydrophobic = analyze_hydrophobic(structure, active_region, surface_region)
				active_hydrophobic.append(normalized_active_region_hydrophobic)
				surface_hydrophobic.append(normalized_surface_region_hydrophobic)
				
				active_fraction, surface_fraction, overall_fraction = analyze_aa(structure, active_region, surface_region)
				
				active_charge_fraction.append(active_fraction["charge"])
				active_polar_fraction.append(active_fraction["polar"])
				active_hydrophobic_fraction.append(active_fraction["hydrophobic"])
				
				surface_charge_fraction.append(surface_fraction["charge"])
				surface_polar_fraction.append(surface_fraction["polar"])
				surface_hydrophobic_fraction.append(surface_fraction["hydrophobic"])
				
				overall_charge_fraction.append(overall_fraction["charge"])
				overall_polar_fraction.append(overall_fraction["polar"])
				overall_hydrophobic_fraction.append(overall_fraction["hydrophobic"])
		
		all_active_area.append(active_hydrophobic)
		all_surface_area.append(surface_hydrophobic)
		
		all_active_charge_fraction.append(active_charge_fraction)
		all_surface_charge_fraction.append(surface_charge_fraction)
		all_overall_charge_fraction.append(overall_charge_fraction)

		
		all_active_polar_fraction.append(active_polar_fraction)
		all_surface_polar_fraction.append(surface_polar_fraction)
		all_overall_polar_fraction.append(overall_polar_fraction)

		all_active_hydrophobic_fraction.append(active_hydrophobic_fraction)
		all_surface_hydrophobic_fraction.append(surface_hydrophobic_fraction)
		all_overall_hydrophobic_fraction.append(overall_hydrophobic_fraction)

	# input()
	from scipy.stats import ranksums
	# active_area
	active_statistic, active_p_value = ranksums(all_active_area[0], all_active_area[1])
	print("active_area")
	print("T 统计值：", active_statistic)
	print("P 值：", active_p_value)
	print("嗜热：", np.mean(all_active_area[0]))
	print("非嗜热：", np.mean(all_active_area[1]))
	write2csv(csv_path=r"D:\subject\active\1-qProtein\data\enzymes\GH11\4_RetAnalyzing\shell\active_area.csv", data=all_active_area)
	print()
	
	# surface_area
	surface_statistic, surface_p_value = ranksums(all_surface_area[0], all_surface_area[1])
	print("surface_area")
	print("T 统计值：", surface_statistic)
	print("P 值：", surface_p_value)
	print("嗜热：", np.mean(all_surface_area[0]))
	print("非嗜热：", np.mean(all_surface_area[1]))
	write2csv(csv_path=r"D:\subject\active\1-qProtein\data\enzymes\GH11\4_RetAnalyzing\shell\surface_area.csv", data=all_surface_area)
	print()
	
	# active_charge_fraction
	active_charge_fraction_statistic, active_charge_fraction_p_value = ranksums(all_active_charge_fraction[0], all_active_charge_fraction[1])
	print("active_charge_fraction")
	print("T 统计值：", active_charge_fraction_statistic)
	print("P 值：", active_charge_fraction_p_value)
	print("嗜热：", np.mean(all_active_charge_fraction[0]))
	print("非嗜热：", np.mean(all_active_charge_fraction[1]))
	write2csv(csv_path=r"D:\subject\active\1-qProtein\data\enzymes\GH11\4_RetAnalyzing\shell\active_charge_fraction.csv", data=all_active_charge_fraction)
	print()
	
	# surface_charge_fraction
	surface_charge_fraction_statistic, surface_charge_fraction_p_value = ranksums(all_surface_charge_fraction[0], all_surface_charge_fraction[1])
	print("surface_charge_fraction")
	print("T 统计值：", surface_charge_fraction_statistic)
	print("P 值：", surface_charge_fraction_p_value)
	print("嗜热：", np.mean(all_surface_charge_fraction[0]))
	print("非嗜热：", np.mean(all_surface_charge_fraction[1]))
	write2csv(csv_path=r"D:\subject\active\1-qProtein\data\enzymes\GH11\4_RetAnalyzing\shell\surface_charge_fraction.csv", data=all_surface_charge_fraction)
	print()
	
	# overall_charge_fraction
	overall_charge_fraction_statistic, overall_charge_fraction_p_value = ranksums(
			all_overall_charge_fraction[0], all_overall_charge_fraction[1]
			)
	print("overall_charge_fraction")
	print("T 统计值：", overall_charge_fraction_statistic)
	print("P 值：", overall_charge_fraction_p_value)
	print("嗜热：", np.mean(all_overall_charge_fraction[0]))
	print("非嗜热：", np.mean(all_overall_charge_fraction[1]))
	write2csv(csv_path=r"D:\subject\active\1-qProtein\data\enzymes\GH11\4_RetAnalyzing\shell\overall_charge_fraction.csv", data=all_overall_charge_fraction)
	print()
	
	# active_polar_fraction
	active_polar_fraction_statistic, active_polar_fraction_p_value = ranksums(
			all_active_polar_fraction[0], all_active_polar_fraction[1]
			)
	print("active_polar_fraction")
	print("T 统计值：", active_polar_fraction_statistic)
	print("P 值：", active_polar_fraction_p_value)
	print("嗜热：", np.mean(all_active_polar_fraction[0]))
	print("非嗜热：", np.mean(all_active_polar_fraction[1]))
	write2csv(csv_path=r"D:\subject\active\1-qProtein\data\enzymes\GH11\4_RetAnalyzing\shell\active_polar_fraction.csv", data=all_active_polar_fraction)
	print()
	
	# surface_polar_fraction
	surface_polar_fraction_statistic, surface_polar_fraction_p_value = ranksums(
			all_surface_polar_fraction[0], all_surface_polar_fraction[1]
			)
	print("surface_polar_fraction")
	print("T 统计值：", surface_polar_fraction_statistic)
	print("P 值：", surface_polar_fraction_p_value)
	print("嗜热：", np.mean(all_surface_polar_fraction[0]))
	print("非嗜热：", np.mean(all_surface_polar_fraction[1]))
	write2csv(csv_path=r"D:\subject\active\1-qProtein\data\enzymes\GH11\4_RetAnalyzing\shell\surface_polar_fraction.csv", data=all_surface_polar_fraction)
	print()
	
	# overall_polar_fraction
	overall_polar_fraction_statistic, overall_polar_fraction_p_value = ranksums(
			all_overall_polar_fraction[0], all_overall_polar_fraction[1]
			)
	print("overall_polar_fraction")
	print("T 统计值：", overall_polar_fraction_statistic)
	print("P 值：", overall_polar_fraction_p_value)
	print("嗜热：", np.mean(all_overall_polar_fraction[0]))
	print("非嗜热：", np.mean(all_overall_polar_fraction[1]))
	write2csv(csv_path=r"D:\subject\active\1-qProtein\data\enzymes\GH11\4_RetAnalyzing\shell\overall_polar_fraction.csv", data=all_overall_polar_fraction)
	print()
	
	# active_hydrophobic_fraction
	active_hydrophobic_fraction_statistic, active_hydrophobic_fraction_p_value = ranksums(
			all_active_hydrophobic_fraction[0], all_active_hydrophobic_fraction[1]
			)
	print("active_hydrophobic_fraction")
	print("T 统计值：", active_hydrophobic_fraction_statistic)
	print("P 值：", active_hydrophobic_fraction_p_value)
	print("嗜热：", np.mean(all_active_hydrophobic_fraction[0]))
	print("非嗜热：", np.mean(all_active_hydrophobic_fraction[1]))
	write2csv(csv_path=r"D:\subject\active\1-qProtein\data\enzymes\GH11\4_RetAnalyzing\shell\active_hydrophobic_fraction.csv", data=all_active_hydrophobic_fraction)
	print()
	
	# surface_hydrophobic_fraction
	surface_hydrophobic_fraction_statistic, surface_hydrophobic_fraction_p_value = ranksums(
			all_surface_hydrophobic_fraction[0], all_surface_hydrophobic_fraction[1]
			)
	print("surface_hydrophobic_fraction")
	print("T 统计值：", surface_hydrophobic_fraction_statistic)
	print("P 值：", surface_hydrophobic_fraction_p_value)
	print("嗜热：", np.mean(all_surface_hydrophobic_fraction[0]))
	print("非嗜热：", np.mean(all_surface_hydrophobic_fraction[1]))
	write2csv(csv_path=r"D:\subject\active\1-qProtein\data\enzymes\GH11\4_RetAnalyzing\shell\surface_hydrophobic_fraction.csv", data=all_surface_hydrophobic_fraction)
	print()
	
	# overall_hydrophobic_fraction
	overall_hydrophobic_fraction_statistic, overall_hydrophobic_fraction_p_value = ranksums(
			all_overall_hydrophobic_fraction[0], all_overall_hydrophobic_fraction[1]
			)
	print("overall_hydrophobic_fraction")
	print("T 统计值：", overall_hydrophobic_fraction_statistic)
	print("P 值：", overall_hydrophobic_fraction_p_value)
	print("嗜热：", np.mean(all_overall_hydrophobic_fraction[0]))
	print("非嗜热：", np.mean(all_overall_hydrophobic_fraction[1]))
	write2csv(csv_path=r"D:\subject\active\1-qProtein\data\enzymes\GH11\4_RetAnalyzing\shell\overall_hydrophobic_fraction.csv", data=all_overall_hydrophobic_fraction)
	print()
	
	# hydrophobic_distribution
	thermo_active = all_active_area[0]
	non_thermo_active = all_active_area[1]
	thermo_active.extend(non_thermo_active)
	
	thermo_surface = all_surface_area[0]
	non_thermo_active = all_surface_area[1]
	thermo_surface.extend(non_thermo_active)
	region_statistic, region_p_value = ranksums(thermo_active, thermo_surface)
	print("hydrophobic distribution")
	print("T 统计值：", region_statistic)
	print("P 值：", region_p_value)
	print("active：", np.mean(thermo_active))
	print("surface：", np.mean(thermo_surface))
	write2csv(csv_path=r"D:\subject\active\1-qProtein\data\enzymes\GH11\4_RetAnalyzing\shell\hydrophobic_distribution.csv", data=[thermo_active, thermo_surface])

template_name = "therm_1TE1_2_Q9HFH0_temperature_55_PDB_1te1_B"
template_shell_file = r"D:\subject\active\1-qProtein\data\enzymes\GH11\2_StrucMapping\therm_1TE1_2_Q9HFH0_temperature_55_PDB_1te1_B_first_shell.txt"
origin_fasta_file = r"D:\subject\active\1-qProtein\data\enzymes\GH11\1_dbBuilder\gh11.fasta"
align_fasta_file = r"D:\subject\active\1-qProtein\data\enzymes\GH11\1_dbBuilder\gh11_usalign.fasta"
structure_dir = [r"D:\subject\active\1-qProtein\data\enzymes\GH11\2_StrucMapping\positive",
                 r"D:\subject\active\1-qProtein\data\enzymes\GH11\2_StrucMapping\negative"]
run(template_name, template_shell_file, origin_fasta_file, align_fasta_file, structure_dir)