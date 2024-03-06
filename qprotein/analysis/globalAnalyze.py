"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/07/30

# Description: Service of characterizing overall structure for a given protein.
# ------------------------------------------------------------------------------
"""
import os

import numpy as np
from qprotein.characterization import disulfide_bond
from qprotein.characterization import hbond
from qprotein.characterization import hydrophobic
from qprotein.characterization import salt_bridge


def statistic_disulfide(atom_array):
    return {"number": len(disulfide_bond.analyze(structure=atom_array))}


def statistic_hbond(structure_path):
    hbond_feature = hbond.analyze(structure_path=structure_path)
    hbond_num = len(hbond_feature)
    avg_hbond_angle = np.mean([i[0] for i in hbond_feature])
    avg_hbond_length = np.mean([i[1] for i in hbond_feature])
    return {"avg_angle": avg_hbond_angle, "avg_length": avg_hbond_length, "number": hbond_num}


def statistic_salt(atom_array):
    salt_feature = salt_bridge.analyze(atom_array=atom_array)
    salt_num = len(salt_feature)
    avg_salt_length = np.mean([i[1] for i in salt_feature])
    return {"avg_salt_length": avg_salt_length, "number": salt_num}


def statistic_hydrophobic(atom_array):
    hydrophobic_feature = hydrophobic.analyze(atom_array=atom_array)
    sum_hydrophobic_area = hydrophobic_feature['sum_area']
    max_hydrophobic_area = hydrophobic_feature['max_area']
    num_hydrophobic = len(hydrophobic_feature['cluster_res'])
    return {"sum_area": sum_hydrophobic_area, "max_area": max_hydrophobic_area, "number": num_hydrophobic}


def analyze(structure_path, atom_array):
    structure_features = {}
    hydrophobic_result = statistic_hydrophobic(atom_array)
    hbond_result = statistic_hbond(structure_path)
    salt_result = statistic_salt(atom_array)
    disulfide_result = statistic_disulfide(atom_array)

    structure_features['hydrophobic_cluster'] = hydrophobic_result
    structure_features['disulfide_bond'] = disulfide_result
    structure_features['hbond'] = hbond_result
    structure_features['salt_bond'] = salt_result
    return structure_features

