"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/05/19

# Description: 
# ------------------------------------------------------------------------------
"""
import os
from Bio import Entrez

work_dir = r"D:\subject\active\1-qProtein\data\enzymes\GH48\1_preprocessing"
# 打开 NCBI Entrez 数据库的访问权限
Entrez.email = "bj600800@gmail.com" # 请写上你自己的电子邮件地址以遵守NCBI的规则
Entrez.api_key = "ec940f6e9e78a1563b00cb2ecc18db028d08" # 如果有NCBI API key请填入，这样可以提高访问速度和配额限制

with open(os.path.join(work_dir, "2_acc_GH48.txt"), "r") as f:
    id_list = [i.rstrip() for i in f.readlines()]
 
# 使用 efetch 函数批量获取序列
handle = Entrez.efetch(db="protein", id=id_list, rettype="fasta", retmode="text")
sequences = handle.read()


with open(os.path.join(work_dir, "3_GH48_sequence.fasta"), "w") as f:
    f.writelines(sequences)
