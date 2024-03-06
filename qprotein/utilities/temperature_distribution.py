"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/10/21

# Description: 
# ------------------------------------------------------------------------------
"""
import os


def get_temperature(dir_path):
	files = os.listdir(dir_path)
	temperature = []
	for i in files:
		t = os.path.splitext(i)[0].split("_")
		print(t)
		if "PDB" in t:
			temperature.append(t[-4])
		elif "." in t[-1]:
			temperature.append(t[-3])
		else:
			temperature.append(t[-1])
	return temperature

dir_path = r"D:\subject\active\1-qProtein\data\enzymes\GH11\2_StrucMapping\positive"
positive_temperature = get_temperature(dir_path)

dir_path2 = r"D:\subject\active\1-qProtein\data\enzymes\GH11\2_StrucMapping\negative"
negative_temperature = get_temperature(dir_path2)
positive_temperature.extend(negative_temperature)
temperature = [str(i)+"\n" for i in positive_temperature]

with open(r"D:\subject\active\1-qProtein\data\enzymes\GH11\temperature.txt", "w") as f:
	f.writelines(temperature)