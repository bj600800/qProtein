"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2024/01/11

# Description: 
# ------------------------------------------------------------------------------
"""
import http.client
import os
from Bio import Entrez
from tqdm import tqdm

def fetch(entry_file):
	with open(entry_file, "r") as f:
		id_list = [i.rstrip().split("\t")[-2] for i in f.readlines()]

	Entrez.email = ""  # 请写上你自己的电子邮件地址以遵守NCBI的规则
	Entrez.api_key = ""  # 如果有NCBI API key请填入，这样可以提高访问速度和配额限制
	
	sequences = []
	batch_size = 50
	total_query = len(id_list)
	for i in range(0, total_query, batch_size):
		try:
			handle = Entrez.efetch(db="protein", id=id_list[i:i+batch_size], rettype="fasta", retmode="text")
			batch_sequences = handle.read()
			sequences.extend(batch_sequences)
		except:
			pass
	

	return sequences

if __name__ == '__main__':
	work_dir = r"D:\subject\active\2-GetThermo\data\TIM"
	tim_dirs = ['GH17', 'GH26', 'GH30', 'GH35',
	            'GH39', 'GH42', 'GH50', 'GH51', 'GH53', 'GH59', 'GH72', 'GH79',
	            'GH86', 'GH113', 'GH128', 'GH140', 'GH147', 'GH148', 'GH157',
	            'GH158', 'GH164', 'GH167', 'GH169', 'GH173']
	for d_name in tqdm(tim_dirs):
		entry_file = os.path.join(work_dir, d_name, d_name+'.txt')
		fasta_file = os.path.join(work_dir, d_name, d_name+'.fasta')
		if not os.path.exists(fasta_file):
			try:
				print()
				print(d_name)
				sequences = fetch(entry_file)
				with open(fasta_file, "w") as f:
					f.writelines(sequences)
			except ValueError as ve:
				print(f"ValueError: {ve}")
			except http.client.IncompleteRead as ire:
				print(f"IncompleteRead: {ire}")
