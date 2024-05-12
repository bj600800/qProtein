"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/09/19

# Description: Predict structure for non-structure sequences from sql.
# ------------------------------------------------------------------------------
"""
import os
import numpy as np
import subprocess

from Bio.PDB import PDBParser
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)

def _get_bfactor(pdb_file):
	parser = PDBParser(QUIET=True)
	structure = parser.get_structure("structure", pdb_file)
	b_factors = [atom.get_bfactor() for atom in structure.get_atoms()]
	return np.mean(b_factors)

def run(fasta_file, structure_dir, esm_script, esm_dir):
	cmd = ["python", esm_script, "--fasta", fasta_file, "--pdb", structure_dir, '-m', esm_dir]
	result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	logger.info(result.stdout.decode())
	for query_name in os.listdir(structure_dir):
		pdb_path = os.path.join(structure_dir, query_name)
		if os.path.exists(pdb_path):
			plddt = _get_bfactor(pdb_path)
			if plddt < 70:
				os.remove(pdb_path)

					

	