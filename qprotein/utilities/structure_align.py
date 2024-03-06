"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/07/22

# Description: align structure with TMalign
# ------------------------------------------------------------------------------
"""
import subprocess

def structure_align(struct1, struct2):
	tm_align_cmd = ["tmalign", "struct1", "struct2"]