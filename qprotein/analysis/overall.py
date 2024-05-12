"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2024/05/10

# Description: 
# ------------------------------------------------------------------------------
"""
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)

def save_feature(feature_file, feature_dict):
	with open(feature_file, "w") as f:
		f.write("Protein_name,Hydrophobic,Hbond,Salt_bridge,Disulfide\n")
		for name, feature in feature_dict.items():
			f.write(name+","+str(feature["hydrophobic"]["sum_area"])+","+str(feature["hbond"])
			        +","+str(feature["saltbridge"])+","+str(feature["disulfide"])+"\n")
	logger.info(f"Overall feature saved: {feature_file}")