"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/04/14

# Description: Step 2. Get structures from AlphaFold structure database.
# ------------------------------------------------------------------------------
"""
import json
import os
import threading
import multiprocessing
from queue import Queue
import time
import random
import tempfile
from tqdm import tqdm
import requests
from requests.adapters import Retry
import numpy as np
from pdbecif.mmcif_io import CifFileWriter
import biotite.structure.io as strucio
import biotite.structure as struc
from biotite.structure.io.pdb import PDBFile

from qprotein.preprocessing.sqlite3_searcher import SqlSearch
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


def start_request_session():
	retries = Retry(total=3, backoff_factor=0.5, status_forcelist=[408, 429, 500, 502, 503, 504])
	session = requests.Session()
	session.mount('http://', requests.adapters.HTTPAdapter(max_retries=retries))
	session.mount('https://', requests.adapters.HTTPAdapter(max_retries=retries))
	return session


class AFDBMapper(SqlSearch):
	def __init__(self, sql_db, folder):
		self.sql_db = sql_db
		self.positive_folder = os.path.join(folder, "positive")
		self.negative_folder = os.path.join(folder, "negative")
		if not os.path.exists(self.positive_folder):
			os.mkdir(self.positive_folder)
		if not os.path.exists(self.negative_folder):
			os.mkdir(self.negative_folder)
			
	def get_uniprot_id(self):
		cursor = self.connect_sql(self.sql_db)
		sql_cmd = "SELECT query_name, accession_id, subject_start, subject_end, label From results_summary " \
				  "WHERE source == 'AFDB'"
		uniprot_ret = self.fetch_results(cursor, sql_cmd)
		return uniprot_ret
	
	def get_exist_query(self):
		exist_structure = [file for folder in [self.positive_folder, self.negative_folder] for file in
		                   os.listdir(folder) if
		                   os.path.isfile(os.path.join(folder, file)) and os.path.getsize(
				                   os.path.join(folder, file)
				                   ) > 0]
		exist_query_name = [os.path.splitext(item)[0] for item in exist_structure]
		return exist_query_name

	def get_AFDB(self, query_info):
		query_name = query_info[0]
		uniprot_id = query_info[1]
		start_residue = query_info[2]
		end_residue = query_info[3]
		label = query_info[4]
		headers = {
			"User-Agent": "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.25 Safari/537.36 Core/1.70.3861.400 QQBrowser/10.7.4313.400",
			"From": "bj600800@gmail.com"  # ALLWAYS TELLs WHO YOU ARE
		}

		api_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
		session = start_request_session()
		try:
			response = session.get(api_url, headers=headers)
			pdb_string = response.text
			if response.status_code == 200:
				with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".pdb") as temp:
					temp.write(pdb_string)
					temp_pdb_filepath = temp.name
					temp.flush()
				
				# 从PDB文件获取结构信息和bfactor信息
				file = PDBFile.read(temp_pdb_filepath)
				bfactor = file.get_b_factor(model=1)
				structure = strucio.load_structure(temp_pdb_filepath)
				structure.set_annotation("b_factor", bfactor)

				# renumber residue id for the wrong residue id
				structure = struc.renumber_res_ids(structure, start=1)
				
				# make path
				if label == "positive":
					clean_file_path = os.path.join(self.positive_folder, f"{query_name}.pdb")
				else:
					clean_file_path = os.path.join(self.negative_folder, f"{query_name}.pdb")
				
				# filter given start and end position residues
				mask_segment = np.isin(structure.res_id, list(range(start_residue, end_residue + 1)))
				structure = structure[mask_segment]

				# renumber residue ids after segmenting structure
				structure = struc.renumber_res_ids(structure, start=1)
				plddt = np.mean(structure.b_factor)
				# save structure
				if plddt > 70:
					strucio.save_structure(clean_file_path, structure)
					os.remove(temp_pdb_filepath)
					afdb_info = [query_name, clean_file_path]
					print("Structure path: ", clean_file_path)
					# return afdb_info
		except:
			raise
	
	def run(self):
		query_info = []
		uniprot_ret = self.get_uniprot_id()
		exist_query_name = self.get_exist_query()
		uniprot_query = [result for result in uniprot_ret if result[0] not in exist_query_name]
		for query in uniprot_query:
			query_info.append(self.get_AFDB(query_info=query))
		# return query_info

if __name__ == '__main__':
	structure_folder = r"D:\subject\active\1-qProtein\data\enzymes\GH5_2\2_StrucMapping"
	summary_sql_path = r"D:\subject\active\1-qProtein\data\enzymes\GH5_2\1_preprocessing\qprotein_results.db"
	test = AFDBMapper(sql_db=summary_sql_path, folder=structure_folder)
	test.run()
