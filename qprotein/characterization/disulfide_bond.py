"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/07/11

# Description: Service of disulfide bond calculation.
The implementation of the following detect_disulfide_bonds function is based on the code by Sphinx-Gallery
Original code source: https://www.biotite-python.org/examples/gallery/structure/disulfide_bonds.html
The modified code is used to detect disulfide bonds in protein structures and returns the frequency of disulfide bonds.
# ------------------------------------------------------------------------------
"""
import numpy as np
import biotite.structure as struc


def analyze(structure, distance=2.05, distance_tol=0.05, dihedral=90, dihedral_tol=15):
    disulfide_bonds = []
    sulfide_mask = (structure.res_name == "CYS") & \
                   (structure.atom_name == "SG")
    cell_list = struc.CellList(
        structure,
        cell_size=distance + distance_tol,
        selection=sulfide_mask
    )
    for sulfide_i in np.where(sulfide_mask)[0]:
        potential_bond_partner_indices = cell_list.get_atoms(
            coord=structure.coord[sulfide_i],
            radius=distance + distance_tol
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
                bond_tuple = tuple(sorted((sulfide_i, sulfide_j)))
                if bond_tuple not in disulfide_bonds:
                    disulfide_bonds.append(bond_tuple)
    return disulfide_bonds

if __name__ == '__main__':
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
    # csv_path = r"D:\subject\active\1-qProtein\data\enzymes\GH8\3_StrucAnalyzing\dsbond.csv"
    #
    # all_ret = []
    # for sdir in structure_dir:
    #     one_ret = []
    #     for struct in tqdm(os.listdir(sdir)):
    #         structure_path = os.path.join(sdir, struct)
    #         structure = strucio.load_structure(structure_path)
    #         seq_len = max(structure.res_id)
    #         out = analyze(structure)
    #         one_ret.append(len(out)/seq_len)
    #     all_ret.append(one_ret)
    # write2csv(csv_path=csv_path, data=all_ret)
    # from scipy.stats import ranksums
    # statistic, p_value = ranksums(all_ret[0], all_ret[1])
    #
    # # 输出结果
    # print("T 统计值：", statistic)
    # print("P 值：", p_value)
    # print(np.mean(all_ret[0]), np.mean(all_ret[1]))
    
    
    dir_path = r"D:\subject\active\3-PETase\data\test"
    import biotite.structure.io as strucio
    import os
    for file in os.listdir(dir_path):
        structure_path = os.path.join(dir_path, file)
        print(file)
        structure = strucio.load_structure(structure_path)
        out = analyze(structure)
        seq_len = max(structure.res_id)
        print(len(out))
        print()