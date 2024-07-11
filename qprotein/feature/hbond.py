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

import biotite.structure as struc
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
    res_id, res_name = struc.get_residues(stack)
    length = len(res_id)
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
    hbond_results = len([tuple(i) for i in zip_results])/length
    # for i in hbond_results:
    #     print('+'.join([str(j) for j in list(stack[i[2]].res_id)]))
    return hbond_results


def analyze(structure_path):
    pdb_addh_temp_file_path = add_h(structure_path)
    hbond_results = calc_hbond(structure_path, pdb_addh_temp_file_path)
    os.remove(pdb_addh_temp_file_path)
    return hbond_results
