"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/04/21

# Description: quantify structures
# ------------------------------------------------------------------------------
"""
import os.path
import sys
import argparse
sys.path.append("..")

from qprotein.structure.StructureAnalyzer import multiprocess
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)

parser = argparse.ArgumentParser(description='Quantify structures')
parser.add_argument('--task_dir', required=True, help='Specific the task directory')
args = parser.parse_args()


def run():
    task_dir = args.task_dir
    struct_dir = os.path.join(task_dir, 'structure')
    output_file = os.path.join(task_dir, 'structure_results.csv')
    multiprocess(struct_dir, output_file)


run()
