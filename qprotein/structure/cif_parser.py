#!/usr/bin/env python           
# -*- coding:utf-8 -*-          
# @Filename:    cif_parser.py      
# @Author:      Eric Dou        
# @Time:        2022/11/11 15:42 

"""
parser cif structure
"""
import os
import pickle
import numpy as np
from matplotlib import pyplot as plt

def get_pae_plddt(model_files):
    np.set_printoptions(threshold=np.inf)
    out = {}

    for i, name in enumerate(model_files):
        d = pickle.load(open(name, 'rb'))
        print(d.keys())
        input()
        out[f'model_{i+1}'] = {'aligned_confidence_probs': d['aligned_confidence_probs']}

    ave = {}
    for model, exp in out.items():
        ave_model = []
        print(exp)
        for i in exp['exp_resol']['logits']:
            ave_model.append(np.average(i))
        ave[model] = ave_model

    plt.subplot(1, 2, 2)
    plt.title("Predicted LDDT per position")
    for model_name, value in ave.items():
        plt.plot(value, label=model_name)
    plt.legend()
    plt.ylim(-10, 20)
    plt.ylabel("Predicted LDDT")
    plt.xlabel("Positions")
    plt.show()
    input()
    return out

model_path = r'E:\pymol\Lib\site-packages\pmg_tk\startup\pyMulstruct\server\src\test_data_visual\wxmy'
model_files = [os.path.join(model_path, i) for i in os.listdir(model_path) if i[-7:] == 'result_model_1_ptm.pkl']
get_pae_plddt(model_files)
