"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/09/19

# Description: Step 3. Predict structure for non-structure sequences from sql.
# ------------------------------------------------------------------------------
"""
import os
import tempfile
import numpy as np
import subprocess

from biotite.structure.io.pdb import PDBFile
from qprotein.preprocessing.sqlite3_searcher import SqlSearch

class ESMMapper(SqlSearch):
	def __init__(self, sql_db, structure_dir, esm_app):
		self.sql_db = sql_db
		self.structure_dir = structure_dir
		self.esm_app = esm_app

	def get_query_fasta(self):
		cursor = self.connect_sql(self.sql_db)
		sql_cmd = "SELECT query_name, sequence From results_summary " \
		          "WHERE struct_path IS NULL"
		fasta_ret = self.fetch_results(cursor, sql_cmd)
		query_name_list = []
		
		# write to temp_file
		with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta") as temp:
			for ret in fasta_ret:
				query_name = ret[0]
				query_name_list.append(query_name)
				sequence = ret[1]
				temp.write(">"+query_name+"\n")
				temp.write(sequence+"\n")
			temp_fasta_path = temp.name
			temp.flush()
		return temp_fasta_path, query_name_list
	
	def ESMfold(self, fasta_file, query_name_list):
		# cmd = ["python", self.esm_app, "-i", fasta_file, "-o", self.structure_dir]
		# subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		# os.remove(fasta_file)
		for query_name in os.listdir(self.structure_dir):
			pdb_path = os.path.join(self.structure_dir, query_name)
			if os.path.exists(pdb_path):
				# 从PDB文件获取结构信息和bfactor信息，小于70删除
				file = PDBFile.read(pdb_path)
				plddt = np.mean(file.get_b_factor(model=1))
				if plddt < 70:
					os.remove(pdb_path)
					
	def run(self):
		temp_fasta_path, query_name_list = self.get_query_fasta()
		print(temp_fasta_path)
		input()
		self.ESMfold(fasta_file='', query_name_list='')

if __name__ == '__main__':
	sql_db = r"D:\subject\active\1-qProtein\data\enzymes\GH48\1_preprocessing\qprotein_results.db"
	struct_dir = r"C:\Users\bj600\Desktop\esm"
	esm = ESMMapper(sql_db=sql_db, structure_dir=struct_dir, esm_app='')
	esm.run()
	