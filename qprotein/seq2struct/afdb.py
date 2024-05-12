"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/04/14

# Description: Step 2. Get structures from AlphaFold structure database.
# ------------------------------------------------------------------------------
"""
import io
import requests
from requests.adapters import Retry
import numpy as np
from Bio.PDB import PDBParser

from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


def _start_request_session():
	retries = Retry(total=3, backoff_factor=0.5, status_forcelist=[408, 429, 500, 502, 503, 504])
	session = requests.Session()
	session.mount('http://', requests.adapters.HTTPAdapter(max_retries=retries))
	session.mount('https://', requests.adapters.HTTPAdapter(max_retries=retries))
	return session


def _get_bfactor(pdb_string):
	pdb_file = io.StringIO(pdb_string)
	parser = PDBParser(QUIET=True)
	structure = parser.get_structure("structure", pdb_file)
	b_factors = [atom.get_bfactor() for atom in structure.get_atoms()]
	return np.mean(b_factors)

def get_struct(uniprot_id):
	headers = {
		"User-Agent": "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.25 Safari/537.36 Core/1.70.3861.400 QQBrowser/10.7.4313.400",
		"From": "bj600800@gmail.com"  # ALLWAYS TELLs WHO YOU ARE
		}
	
	api_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
	session = _start_request_session()
	response = session.get(api_url, headers=headers)
	if response.status_code == 200:
		pdb_string = response.text
		# 从PDB文件获取结构信息和bfactor信息
		plddt = _get_bfactor(pdb_string)
		if plddt > 70:
			return pdb_string
	else:
		logger.info(f"Seq2Struct uniprot ID failed: {uniprot_id}")



# if __name__ == '__main__':
# 	structure_folder = r"D:\subject\active\1-qProtein\code\test\structure"
# 	get_struct("A0A067XRW0", structure_folder)
