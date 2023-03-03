#!/usr/bin/env python           
# -*- coding:utf-8 -*-          
# @Filename:    pkl_file.py      
# @Author:      Eric Dou        
# @Time:        2022/5/5 9:04 

"""
Reassign plddt value for every residue.

Every structure was valued through a single float number predicted by AlphaFold 2.0.0. But the plddts
could be reassigned to every residue from a protein structure, which enable users to analysis more confidently.

Here is the code logic:
    1. Locate output files of every target.
    2. Get plddt values for every structure from corresponding pkl file.
    3. Replace B factor for every pdb structure file.
"""

import os
import json
import pickle
from Bio.PDB import PDBParser, PDBIO
import warnings
import shutil
import plotly.express as px

warnings.filterwarnings("ignore")
io = PDBIO()
p = PDBParser()


def reass_plddt(target):
    """
    reassign plddt value to every structure in target directory.
    @param target:
    @return:
    """
    # reassign plddt for every atom of every structure
    for num in range(1, 6):
        with open(os.path.join(target, f'result_model_{num}.pkl'), 'rb') as pkl:
            data = pickle.load(pkl)
            plddt = data["plddt"]
            structure = p.get_structure(f'relaxed_model_{num}', os.path.join(target, f'relaxed_model_{num}.pdb'))
            atoms = structure.get_atoms()
            for atom in atoms:
                # atom.set_bfactor(plddt[atom.get_full_id()[3][1] - 1]) # Original value
                atom.set_bfactor(100 - plddt[atom.get_full_id()[3][1] - 1]) # reverse value
            io.set_structure(structure)
            io.save(os.path.join(target, f'relaxed_model_bf_{num}.pdb'), preserve_atom_numbering=True)

    # rename relaxed_model_bf_{num}.pdb to ranked_bf_{idx}.pdb
    with open(os.path.join(target, r'ranking_debug.json'), 'r') as j:
        ranking = json.load(j)
        order = ranking["order"]
        for idx, model_num in enumerate(order):
            shutil.copyfile(os.path.join(target, f'relaxed_model_bf_{order[idx][-1]}.pdb'),
                            os.path.join(target, f'ranked_bf_{idx}.pdb'))


def main(wd_path):
    """
    iterate all targets.
    @param wd_path:
    @return:
    """
    target = [i for i in os.listdir(wd_path) if os.path.isdir(os.path.join(wd_path, i))
              if os.path.isfile(os.path.join(wd_path, i, "ranked_0.pdb"))]
    length = len(target)
    for id, i in enumerate(target):
        reass_plddt(os.path.join(wd_path, i))
        print('Progress:\t',end='')
        print((id+1)/length)

    df = px.data.tips()
    fig = px.scatter(df, x="total_bill", y="tip", color="size",
                     title="Numeric 'size' values mean continuous color")

    fig.show()
    # reass_plddt(r'C:\Users\user\Desktop\NP_001032581.1') # test


# 1. Test the script(af_patch.py) normal or abnormal. 2. In case of duplicated call.
# Explanation:
# Once the script was imported by other codes, the value of __name__ will be replaced by the script's name(af_patch),
# instead of __main__.

if __name__ == "__main__":
    # For only one data directory path
    wd_path = r'./zh_out'
    main(wd_path)

    # # For multiple directory paths
    # for i in range(1,5):
    #     wd_path = f'E:\\af2\data\gcf_{i}\\task'
    #     main(wd_path)

