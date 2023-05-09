"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/04/15

# Description: move protein structures for paper.
# ------------------------------------------------------------------------------
"""
import os
import shutil
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(description='Copy structure files')
parser.add_argument('--task_name', required=True, help='Task name')
parser.add_argument('--task_dir', required=True, help='Specific the task directory')
args = parser.parse_args()

task_name = args.task_name
root_path = args.task_dir
work_path = os.path.join(root_path, task_name)
structure_path = os.path.join(work_path, 'structure')
file = os.listdir(structure_path)
high_plddt_dir = os.path.join(work_path, 'high_structure')
good_plddt_dir = os.path.join(work_path, 'confident_structure')
if not os.path.exists(high_plddt_dir):
    os.mkdir(high_plddt_dir)
if not os.path.exists(good_plddt_dir):
    os.mkdir(good_plddt_dir)

for i in tqdm(file):
    plddt = os.path.splitext(i)[0].split("#")[-1]
    if float(plddt) >= 90:
        source = os.path.join(structure_path, i)
        target = os.path.join(high_plddt_dir, i)
        shutil.copy(source, target)
    else:
        source = os.path.join(structure_path, i)
        target = os.path.join(good_plddt_dir, i)
        shutil.copy(source, target)

