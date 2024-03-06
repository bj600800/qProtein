"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2024/03/06

# Description: Make sequence databases for blast.
# ------------------------------------------------------------------------------
"""

import os
from multiprocessing import cpu_count

from diamond4py import Diamond

from qprotein.utilities import logger
logger = logger.setup_log(name=__name__)

# Count cpu threads
cpu_count = cpu_count()

def check_file_exist(file):
	# Check diamond database existence
	return os.path.isfile(file)

def makedb(db_dir, pdb_fasta, sprot_fasta, trembl_fasta):
	fasta_list = [pdb_fasta, sprot_fasta, trembl_fasta]
	dmnd_list = [os.path.join(db_dir, 'pdb.dmnd'),
	             os.path.join(db_dir, 'sprot.dmnd'),
	             os.path.join(db_dir, 'trembl.dmnd')]
	
	zip_list = zip(dmnd_list, fasta_list)
	for dmnd, fasta in zip_list:
		if not check_file_exist(dmnd):
			diamond = Diamond(
					database=dmnd,
					n_threads=cpu_count
					)
			diamond.makedb(fasta)
		else:
			logger.info(f"Database already exists, skipping the database building process for {os.path.basename(dmnd)}")

	
if __name__ == '__main__':
    db_dir = r"D:\subject\active\1-qProtein\code\db"
    pdb_fasta = r"D:\subject\active\1-qProtein\code\db\pdb.fasta"
    sprot_fasta = r"D:\subject\active\1-qProtein\code\db\sprot.fasta"
    trembl_fasta = r"D:\subject\active\1-qProtein\code\db\trembl.fasta"
    makedb(db_dir, pdb_fasta, sprot_fasta, trembl_fasta)
	