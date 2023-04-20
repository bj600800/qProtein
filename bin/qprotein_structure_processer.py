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

parser.add_argument('--work_dir', required=True, help='All task dirs should be here')
parser.add_argument('--task_name', required=True, help='Name the task')

args = parser.parse_args()


def crawl_structure(work_dir, task_name):
    task_dir = os.path.join(work_dir, task_name)
    structure_dir_path = os.path.join(task_dir, 'structure')
    sql_db = os.path.join(task_dir, 'qprotein_results.db')
    log_file = os.path.join(task_dir, '404NotFoundURL.txt')
    map_multiprocess(sql_db, structure_dir_path, log_file)


def run():
    work_dir = args.work_dir
    task_name = args.task_name
    crawl_structure(work_dir, task_name)


run()
