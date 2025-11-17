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
import biotite.structure.io as strucio

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
	length = len(set(atom_array.res_id))
	return len(salt_bridges)/length

