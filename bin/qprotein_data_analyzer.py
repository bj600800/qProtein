"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/05/03

# Description: The service code for qProtein results analysis
# ------------------------------------------------------------------------------
"""
import os
import csv
from qprotein.utilities.Analyze import SearchParam


def search_params(task_dir):
    ident = ["100"]  # ["70", "80", "90", "100"]
    cover = ["100"]  # ["90", "100"]
    for i in ident:
        for c in cover:
            if i == "95" and c == "90":
                return
            param = i + "_200_" + c
            print(param)
            searcher = SearchParam(param=param, task_dir=task_dir)
            structure_attr = searcher.get_structure_attr()
            with open(os.path.join(task_dir, "structure_results_"+param+".csv"), "w", newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['query_name', 'helix', 'sheet', 'loop', 'turn', 'bend', 'bridge', 'hbond_density',
                                 'hbond_avg_energy', 'apolar', 'polar', 'positive', 'negative',
                                 'polar_area', 'apolar_area']
                                )
                writer.writerows(structure_attr)


task_dir = r"D:\subject\active\1-qProtein\data\tibet"
search_params(task_dir)
