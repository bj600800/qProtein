"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/07/19

# Description: Service for hbond analysis.
Pdb2Pqr for add hydrogen atom, biotite_hbond for calculate h-bonds
# ------------------------------------------------------------------------------
"""
import os
import subprocess
import tempfile
import numpy as np

import biotite.structure.io as strucio
from biotite.structure import hbond, angle, distance
import warnings
warnings.filterwarnings("ignore")


def add_h(input_pdb):
    pqr_temp_file = tempfile.NamedTemporaryFile(delete=False)
    pqr_temp_file_path = pqr_temp_file.name

    pdb_addh_temp_file = tempfile.NamedTemporaryFile(delete=False)
    pdb_addh_temp_file_path = pdb_addh_temp_file.name + ".pdb"

    pdb2pqr_cmd = ["pdb2pqr", "--ff=AMBER", "--log-level", "CRITICAL", "--titration-state-method", "propka",
                   "--pdb-output", pdb_addh_temp_file_path, input_pdb, pqr_temp_file_path]
    pqr_temp_file.close()
    os.remove(pqr_temp_file_path)
    subprocess.run(pdb2pqr_cmd, check=True)
    return pdb_addh_temp_file_path


def calc_hbond(pdb, pdb_h):
    stack = strucio.load_structure(pdb)
    stack_h = strucio.load_structure(pdb_h)
    triplets_h = hbond(stack_h)
    donor_acceptor = [(list(stack_h[bond].res_id[::2]), list(stack_h[bond].atom_name[::2])) for bond in triplets_h]
    triplets = []
    for resID_atmN in donor_acceptor:
        res_id = resID_atmN[0]
        atom_name = resID_atmN[1]
        pair_idx = []
        for i in range(0, 2):
            hbond_mask = (stack.res_id == int(res_id[i])) & (stack.atom_name == atom_name[i])
            pair_idx.append(np.where(hbond_mask)[0])
        pair_idx = np.concatenate(pair_idx)
        triplets.append(pair_idx)
    triplets = np.array(triplets)
    donor = stack_h[triplets_h[:, 0]]
    donor_h = stack_h[triplets_h[:, 1]]
    acceptor = stack_h[triplets_h[:, 2]]
    theta = angle(donor, donor_h, acceptor)
    degree = np.rad2deg(theta)
    dist = distance(donor_h, acceptor)
    zip_results = [*zip(degree, dist, triplets)]
    hbond_results = len([tuple(i) for i in zip_results])
    # for i in hbond_results:
    #     print('+'.join([str(j) for j in list(stack[i[2]].res_id)]))
    return hbond_results


def analyze(structure_path):
    pdb_addh_temp_file_path = add_h(structure_path)
    hbond_results = calc_hbond(structure_path, pdb_addh_temp_file_path)
    os.remove(pdb_addh_temp_file_path)
    return hbond_results


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
    def write2txt(txt_path, data):
        with open(txt_path, "w") as f:
            for i in data:
                f.write(i[0]+"\t"+str(i[1])+"\n")
    import os
    from tqdm import tqdm

    # structure_dir = [r"D:\subject\active\1-qProtein\data\enzymes\GH10\2_StrucMapping\positive",
    #                  r"D:\subject\active\1-qProtein\data\enzymes\GH10\2_StrucMapping\negative"]
    # csv_path = r"D:\subject\active\1-qProtein\data\enzymes\GH8\3_StrucAnalyzing\hbond_num.csv"
    #
    # all_ret = []
    # for sdir in structure_dir:
    #     one_ret = []
    #     for struct in tqdm(os.listdir(sdir)):
    #         structure_path = os.path.join(sdir, struct)
    #         print(structure_path)
    #         structure = strucio.load_structure(structure_path)
    #         seq_len = max(structure.res_id)
    #         hbond_feature = analyze(structure_path=structure_path)
    #         input()
    #         hbond_num = len(hbond_feature)
    #         one_ret.append(hbond_num/seq_len)
    structure_dir = r"C:\Users\bj600\Desktop\ms"
    txt_path = r"C:\Users\bj600\Desktop\hbond_num.txt"
    ret = []
    for structure_file in os.listdir(structure_dir):
        structure_path = os.path.join(structure_dir, structure_file)
        structure = strucio.load_structure(structure_path)
        hbond_feature = analyze(structure_path=structure_path)
        hbond_num = len(hbond_feature)
        ret.append([structure_file.split(".")[0], hbond_num])
    write2txt(txt_path=txt_path, data=ret)
    #     all_ret.append(one_ret)
    # write2csv(csv_path=csv_path, data=all_ret)
    #
    # from scipy.stats import ranksums
    #
    # positive = all_ret[0]
    # negative = all_ret[1]
    # statistic, p_value = ranksums(positive, negative)
    #
    # # 输出结果
    # print("T 统计值：", statistic)
    # print("P 值：", p_value)
    # print(np.mean(positive), np.mean(negative))
    # print()
    
    # dir_path = r"D:\subject\active\3-PETase\data\test"
    # for file in os.listdir(dir_path):
    #     structure_path = os.path.join(dir_path, file)
    #     print(file)
    #     structure = strucio.load_structure(structure_path)
    #     out = analyze(structure_path)
    #     seq_len = max(structure.res_id)
    #     # print(len(out)/seq_len)
    #     print(len(out)/seq_len)
    #     print()
