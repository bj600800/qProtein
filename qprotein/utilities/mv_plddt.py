"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/04/15

# Description: 
# ------------------------------------------------------------------------------
"""
import os
import shutil
from tqdm import tqdm

task_name = "manure"
root_path = r"D:\subject\active\1-qProtein\data"
work_path = os.path.join(root_path, task_name)
structure_path = os.path.join(work_path, 'structure')
statistics_file = os.path.join(work_path, 'statistics_plddt.txt')
file = os.listdir(structure_path)
high_plddt_dir = os.path.join(work_path, 'high_plddt')
good_plddt_dir = os.path.join(work_path, 'good_plddt')
if not os.path.exists(high_plddt_dir):
    os.mkdir(high_plddt_dir)
if not os.path.exists(good_plddt_dir):
    os.mkdir(good_plddt_dir)

for i in tqdm(file):
    plddt = os.path.splitext(i)[0].split("#")[2]
    if float(plddt) >= 90:
        source = os.path.join(structure_path, i)
        target = os.path.join(high_plddt_dir, i)
        shutil.copy(source, target)
    else:
        source = os.path.join(structure_path, i)
        target = os.path.join(good_plddt_dir, i)
        shutil.copy(source, target)

