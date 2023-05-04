"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: Perform structure quantification and sql insert
# ------------------------------------------------------------------------------
"""
import os
import sys
sys.path.append("../")

import argparse
from qprotein.structure.StructureMapper import map_multiprocess


parser = argparse.ArgumentParser(description='Crawl query structure with multiprocessing.')

parser.add_argument('--task_dir', required=True, help='Specific the task directory')

args = parser.parse_args()


def run():
    task_dir = args.task_dir
    map_multiprocess(task_dir)


run()
