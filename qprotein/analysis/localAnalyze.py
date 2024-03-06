"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/07/28

# Description: Service of local structure extraction for a given protein.
Investigate the interactions of the local protein structure within the designated spatial zone.
# ------------------------------------------------------------------------------
"""


import os
import numpy as np
import csv
from tqdm import tqdm

from collections import defaultdict
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as pdb

from qprotein.structure.analyzer.methods import disulfide_bond
from qprotein.structure.analyzer.methods import salt_bridge
from qprotein.structure.analyzer.methods import hbond
from qprotein.utilities import align_map


def construct_cube(atomArray, cell_size=5, k_length=1):
    """
    cube: each position of a given atom alpha carbon around 10 A radii. The position is after-aligned.
    :return: cube_list contained atom ids.
    """
    alpha_carbon_mask = atomArray.atom_name == "CA"
    cell_list = struc.CellList(
        atomArray,
        cell_size=cell_size
    )
    k_length = k_length
    carbon_list = np.where(alpha_carbon_mask)[0]
    cube_list = []
    for carbon_idx in range(0, len(carbon_list), k_length):
        target_idx = carbon_list[carbon_idx:carbon_idx+k_length]

        atom_info = np.array([atomArray[atom_idx].coord for atom_idx in target_idx])
        ###
        indices = cell_list.get_atoms_in_cells(atom_info, cell_radius=1)
        # atom id +1
        cube_atom_idx = set([atom for cube in indices for atom in cube if atom+1 != 0])
        cube_list.append(tuple((target_idx, cube_atom_idx)))
    return cube_list


def array_midpoint(arr):
    # 将数组转换为 NumPy 数组
    arr = np.array(arr)
    # 计算每个元素的中位数
    midpoint = np.median(arr, axis=0)
    return midpoint


def deduplicate_bond(cube_list):
    """
    every bond only contained in one cube.
    """
    cube_dict = defaultdict(list)
    # temp obj
    bond_dist_dict = defaultdict(list)

    # construct bond_dist_dict
    for item in cube_list:
        if item[1]:
            key = tuple(item[2])
            bond_dist_dict[key].append([item[0], item[1]])
    for bond, dist in bond_dist_dict.items():
        sorted_data = sorted(dist, key=lambda x: x[1])
        cube_dict[sorted_data[0][0]].append(bond)
    cube_dict = {k: v for k, v in cube_dict.items()}
    return cube_dict


def assign_interaction_with_cube(features_to_distribute, atom_array):
    """bond numbers in every cube"""
    cube_list = construct_cube(atom_array)
    distributed = []
    for i, cube in enumerate(cube_list):
        _cube_feature = []
        for feature_key, feature_val in features_to_distribute.items():
            if feature_key == "salt_bond":
                result_array = np.array([num_pair for num_pair in feature_val if all(num in cube[1] for num in num_pair)])
                if result_array.any():
                    for ret in result_array:
                        alphac = atom_array[cube[0]].coord
                        midpoint = array_midpoint(atom_array[ret].coord)
                        dist = np.linalg.norm(midpoint - alphac)
                        if list(ret) not in distributed:
                            distributed.append((i+1, dist, list(ret)))
            if feature_key == "hbond":
                result_array = [num_pair for num_pair in feature_val if all(num in cube[1] for num in num_pair)]
                if result_array:
                    for ret in result_array:
                        ret = np.array(ret)
                        alphac = atom_array[cube[0]].coord
                        midpoint = array_midpoint(atom_array[ret].coord)
                        dist = np.linalg.norm(midpoint - alphac)
                        if list(ret) not in distributed:
                            distributed.append((i+1, dist, list(ret)))
            if feature_key == "disulfide_bond":
                result_array = np.array(
                        [num_pair for num_pair in feature_val if all(num in cube[1] for num in num_pair)]
                        )
                if result_array.any():
                    for ret in result_array:
                        alphac = atom_array[cube[0]].coord
                        midpoint = array_midpoint(atom_array[ret].coord)
                        dist = np.linalg.norm(midpoint - alphac)
                        if list(ret) not in distributed:
                            distributed.append((i + 1, dist, list(ret)))
    
    distributed = deduplicate_bond(distributed)
    return distributed


def analyze():
    work_dir = r"D:\subject\active\1-qProtein\data\enzymes\cazy_endo-1_4-beta-xylanase\classification"
    struct_dir = [os.path.join(work_dir, "mesophilic_pdb"),
                  os.path.join(work_dir, "thermophilic_pdb")]
    res_id_map = pairwise.run()
    inverted_dict = {outer_key: {inner_value: inner_key for inner_key, inner_value in inner_dict.items()}
                     for outer_key, inner_dict in res_id_map.items()}
    for strucdir in struct_dir[:1]:
        array_dir = []
        for i in os.listdir(strucdir):
            struct_path = os.path.join(strucdir, i)
            structure_name = os.path.splitext(os.path.basename(struct_path))[0]
            atom_array = read_structure(struct_path=struct_path)
            # features_to_distribute = {"salt_bond": salt_bridge.analyze(atom_array),
            #     "hbond": [i[2] for i in hbond.analyze(struct_path)],
            #     "disulfide_bond": disulfide_bond.analyze(atom_array)}
            # features_to_distribute = {"disulfide_bond": disulfide_bond.analyze(atom_array)}
            features_to_distribute = {"hbond": [i[2] for i in hbond.analyze(struct_path)]}
            # features_to_distribute = {"salt_bond": salt_bridge.analyze(atom_array)}
            distributed_bond = assign_interaction_with_cube(features_to_distribute=features_to_distribute, atom_array=atom_array)
            distributed_bond = {res_id_map[structure_name][k]: v for k, v in distributed_bond.items()}
            if distributed_bond:
                for key, values in distributed_bond.items():
                    if key == 73:
                        print(struct_path)
                        print(structure_name)
                        origin_res = inverted_dict[structure_name][key]
                        print("cube num:", key)
                        print("original residue:", origin_res)
                        for v in values:
                            resi_str = '+'.join([str(atom_array[i].res_id) for i in v])
                            cmd = f"select res {resi_str}"
                            print(cmd)
                        print()
            distributed_num = {k: len(v) for k, v in distributed_bond.items()}
            bond_array = [distributed_num.get(i, 0) for i in range(1, 364)]
            array_dir.append(bond_array)
    
        # array = np.array(array_dir)
        # np.savetxt(fname=strucdir+"_hbond_all.csv", X=array, delimiter=',') # construct input data
        # mean_bond = np.mean(array, axis=0)
        # std_bond = np.std(array, axis=0)
        # with open(strucdir+"_hbond.csv", 'w', newline='') as csvfile:
        #     writer = csv.writer(csvfile)
        #     # 写入表头
        #     writer.writerow(['mean', 'std'])
        #     for mean, std in zip(mean_bond, std_bond):
        #         writer.writerow([mean, std])
            
analyze()
