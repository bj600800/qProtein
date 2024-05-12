"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/07/22

# Description: Service of salt-bridge search in protein structure.
# ------------------------------------------------------------------------------
"""
import numpy as np
import biotite.structure as struc


def detect_disulfide_bonds(structure, distance=4):
	pos_label = {
		"HIS": ("ND1", "NE2"),
		"HSD": ("ND1", "NE2"),
		"HSE": ("ND1", "NE2"),
		"HSP": ("ND1", "NE2"),
		"HIE": ("ND1", "NE2"),
		"HIP": ("ND1", "NE2"),
		"HID": ("ND1", "NE2"),
		"LYS": ("NZ",),
		"ARG": ("NH1", "NH2")
		}
	
	neg_label = {
		"ASP": ("OD1", "OD2"),
		"GLU": ("OE1", "OE2")
		}
	salt_bridges = []
	try:
		for pos_res_name, pos_atom_names in pos_label.items():
			for one_pos_atom_name in pos_atom_names:
				pos_mask = (structure.res_name == pos_res_name) & (structure.atom_name == one_pos_atom_name)
				for neg_res_name, neg_atom_names in neg_label.items():
					for one_neg_atom_name in neg_atom_names:
						neg_mask = (structure.res_name == neg_res_name) & (structure.atom_name == one_neg_atom_name)
						for pos_idx in np.where(pos_mask)[0]:
							for neg_idx in np.where(neg_mask)[0]:
								pos_atom = structure[pos_idx]
								neg_atom = structure[neg_idx]
								bond_dist = struc.distance(pos_atom, neg_atom)
								if bond_dist < distance:
									bridge = sorted([pos_idx, neg_idx])
									if bridge not in [i[0] for i in salt_bridges]:
										salt_bridges.append([bridge, bond_dist])
		return salt_bridges
	except:
		raise


def analyze(atom_array):
	salt_bridges = detect_disulfide_bonds(structure=atom_array)
	return len(salt_bridges)

def write2txt(txt_path, data):
    with open(txt_path, "w") as f:
        for i in data:
            f.write(i[0]+"\t"+str(i[1])+"\n")
		    
if __name__ == '__main__':
	import os
	import biotite.structure.io as strucio
	
	structure_dir = r"C:\Users\bj600\Desktop\ms"
	txt_path = r"C:\Users\bj600\Desktop\saltbridge_num.txt"
	ret = []
	for structure_file in os.listdir(structure_dir):
		structure_path = os.path.join(structure_dir, structure_file)
		structure = strucio.load_structure(structure_path)
		salt_feature = analyze(atom_array=structure)
		if salt_feature: ## ???
			salt_num = len(salt_feature)
			ret.append([structure_file.split(".")[0],salt_num])
			write2txt(txt_path=txt_path, data=ret)
	# import csv
	# from itertools import zip_longest
	# def write2csv(csv_path, data):
	#     max_length = max(len(lst) for lst in data)
	#     transposed_data = list(zip_longest(*data, fillvalue=None))
	#
	#     with open(csv_path, 'w', newline='') as csvfile:
	#         writer = csv.writer(csvfile)
	#         # 写入数据
	#         for i, row in enumerate(transposed_data, start=1):
	#             writer.writerow([i] + list(row[:max_length]))
	#
	# import os
	# from tqdm import tqdm
	# import biotite.structure.io as strucio
	#
	# structure_dir = [r"D:\subject\active\1-qProtein\data\enzymes\GH8\2_StrucMapping\positive",
	#                  r"D:\subject\active\1-qProtein\data\enzymes\GH8\2_StrucMapping\negative"]
	# csv_path = r"D:\subject\active\1-qProtein\data\enzymes\GH8\3_StrucAnalyzing\salt_bridge.csv"
	#
	# all_ret = []
	# for sdir in structure_dir:
	#     one_ret = []
	#     for struct in tqdm(os.listdir(sdir)):
	#         structure_path = os.path.join(sdir, struct)
	#         structure = strucio.load_structure(structure_path)
	#         seq_len = max(structure.res_id)
	#         salt_feature = analyze(atom_array=structure)
	#         if salt_feature: ## ???
	# 	        salt_num = len(salt_feature)
	# 	        one_ret.append(salt_num/seq_len)
	#     all_ret.append(one_ret)
	# write2csv(csv_path=csv_path, data=all_ret)
	# from scipy.stats import ranksums
	# statistic, p_value = ranksums(all_ret[0], all_ret[1])
	#
	# # 输出结果
	# print("T 统计值：", statistic)
	# print("P 值：", p_value)
	# print(np.mean(all_ret[0]), np.mean(all_ret[1]))
	# import os
	# import biotite.structure.io as strucio
	# dir_path = r"D:\subject\active\3-PETase\data\test"
	# for file in os.listdir(dir_path):
	# 	structure_path = os.path.join(dir_path, file)
	# 	print(file)
	# 	structure = strucio.load_structure(structure_path)
	# 	out = analyze(structure)
	# 	seq_len = max(structure.res_id)
	# 	print(len(out))
	# 	print()