"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/10/31

# Description: 
# ------------------------------------------------------------------------------
"""
from biotite.structure.io.pdb import PDBFile
import numpy as np
import os


plddt_90 = []
plddt_70 = []
other_plddt = []
dir_path = [r"D:\subject\active\1-qProtein\data\enzymes\GH48\2_StrucMapping\positive",
            r"D:\subject\active\1-qProtein\data\enzymes\GH48\2_StrucMapping\negative"]
for path in dir_path:
	for file in os.listdir(path):
		structure_path = os.path.join(path, file)
		file_obj = PDBFile.read(structure_path)
		plddt = np.mean(file_obj.get_b_factor(model=1))
		if plddt >= 90:
			plddt_90.append((file, plddt))
		elif 70 <= plddt < 90:
			plddt_70.append((file, plddt))
		elif 0 < plddt <= 70:
			other_plddt.append((file, plddt))
			
print(len(other_plddt+plddt_90+plddt_70))
print(len(plddt_90))
print(len(plddt_70))
for i in other_plddt:
	print(i)
