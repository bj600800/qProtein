"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2024/05/10

# Description: Main script of qProtein
# ------------------------------------------------------------------------------
"""
import os.path
import argparse
import subprocess
import time

from tqdm import tqdm
import biotite.structure.io as strucio
from qprotein.seq2struct.sequence import get_sequence, get_id
from qprotein.seq2struct.afdb import get_struct
from qprotein.seq2struct import esmfold
from qprotein.feature import hydrophobic, hbond, salt_bridge, disulfide_bond
from qprotein.analysis import overall, local
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)

#### ARGUMENTS PARSER ####
parser = argparse.ArgumentParser(description='qProtein for structures analysis')

parser.add_argument('--fasta', required=False, help='Input fasta sequences')
parser.add_argument('--id', required=False, help='Input uniprot ID list')

parser.add_argument('--work_dir', required=True, help='Working directory')
parser.add_argument('--pre_pdb', required=False, help='Directory of prepared pdbs')

parser.add_argument('--local', action='store_true', required=False, default=False, help='Analysis mode=local')
parser.add_argument('--template_name', required=False, help='Template model for local analysis')
parser.add_argument('--template_active_res', required=False, help='Active residues of template model, i.e. 33,35,37,64')
parser.add_argument('--dist1', required=False, help='Active region distance')
parser.add_argument('--dist2', required=False, help='Intermediate region distance')

args = parser.parse_args()

#### USER CONFIGURATION ####
# Binary executable of US-align tool
usalign_binary_path = r""
# assert os.path.exists(usalign_binary_path)
# ESMFold prediction script. Needed if the --fasta is chosen.
esm_script = r"/opt/app/esm-main/scripts/fold.py"
# ESMFold directory, looking for checkpoints directory. Needed if the --fasta is chosen.
esm_dir = r"/opt/app/esm-main/"
#### END OF USER CONFIGURATION ####

if args.fasta:
	logger.info(f"Predicting structures by ESMFold")
elif args.id:
	logger.info(f"Crawl structures from AlphaFold structure database")
	
if args.local:
	if not (args.template_name and args.template_active_res and args.dist1 and args.dist2):
		parser.error("When using local analysis mode, you must specify "
					 "--template_name, --template_active_res, --dist1, --dist2")
		parser.print_help()
		exit(1)
#### END OF ARGUMENTS PARSER ####


def sequences(fasta_file):
	return get_sequence(fasta_file)

def uniprot_ids(id_file):
	return get_id(id_file)

def crawl_struct(id_list, structure_folder):
	def search_exist_struct(structure_folder):
		exist_structure = [file for file in os.listdir(structure_folder)
		                   if os.path.isfile(os.path.join(structure_folder, file))
		                   and os.path.getsize(os.path.join(structure_folder, file)) > 0]
		return [os.path.splitext(item)[0] for item in exist_structure]
	
	id_list = list(set(id_list) - set(search_exist_struct(structure_folder)))
	if id_list:
		for uniprot_id in tqdm(id_list):
			pdb_string = get_struct(uniprot_id)
			if pdb_string:
				save_path = os.path.join(structure_folder, uniprot_id + ".pdb")
				with open(save_path, "w") as f:
					f.write(pdb_string)
			
	logger.info(f"Crawled structures: {len(os.listdir(structure_folder))}")

def calc_hydrophobic(pdb_file):
	structure = strucio.load_structure(pdb_file)
	return hydrophobic.analyze(atom_array=structure)

def calc_hbond(pdb_file):
	return hbond.analyze(structure_path=pdb_file)

def calc_saltbridge(pdb_file):
	structure = strucio.load_structure(pdb_file)
	return salt_bridge.analyze(atom_array=structure)

def calc_disulfide(pdb_file):
	structure = strucio.load_structure(pdb_file)
	return disulfide_bond.analyze(structure)
	

def calc_feature(pdb_list):
	logger.info(f"Calculating features, please wait...")

	feature_dict = {}
	for pdb_file in tqdm(pdb_list):
		feature = {}
		pdb_name = os.path.basename(pdb_file).split(".")[0]
		feature["hydrophobic"] = calc_hydrophobic(pdb_file)
		feature["hbond"] = calc_hbond(pdb_file)
		feature["saltbridge"] = calc_saltbridge(pdb_file)
		feature["disulfide"] = calc_disulfide(pdb_file)
		feature_dict[pdb_name] = feature

	return feature_dict


def align_structure(structure_folder):
	filenames = os.listdir(structure_folder)
	name_txt = os.path.join(args.work_dir, 'name.txt')
	with open(name_txt, 'w') as nf:
		for file in filenames:
			nf.write(file+'\n')
	align_output_file = os.path.join(os.path.dirname(structure_folder), "usalign_out.fasta")

	usalign_cmd = [usalign_binary_path, "-dir", structure_folder, name_txt, "-suffix", ".pdb",
	               "-mm", "4"]
	result = subprocess.run(usalign_cmd, check=True, stdout=subprocess.PIPE)
	with open(align_output_file, 'w', newline='') as f:
		f.writelines(result.stdout.decode())
	os.remove(name_txt)
	return align_output_file

	
def main():
	t1 = time.time()
	# wd = r"D:\subject\active\1-qProtein\code\test"
	wd = args.work_dir
	structure_folder = os.path.join(wd, "structure")
	if not os.path.exists(structure_folder):
		os.mkdir(structure_folder)
	# structure_folder = r"D:\subject\active\1-qProtein\code\test\structure"

	if args.fasta:
		# ESMFOLD prediction
		fasta_file = args.fasta
		esmfold.run(fasta_file=fasta_file, structure_dir=structure_folder, esm_script=esm_script, esm_dir=esm_dir)
	
	if args.id:
		# AFDB crawl
		id_file = args.id
		id_list = get_id(id_file)
		crawl_struct(id_list, structure_folder)

	# Using user prepared structures.
	if args.pre_pdb:
		structure_folder = args.pre_pdb

	# Calculate feature
	pdb_list = [os.path.join(structure_folder, i) for i in os.listdir(structure_folder)]
	feature_dict = calc_feature(pdb_list)
	
	
	# Analyze overall feature
	overall_feature_file = os.path.join(wd, "overall_feature.csv")
	overall.save_feature(overall_feature_file, feature_dict)

	if args.local:
		# Analyze local feature
		hydrophobic_feature = {protein_id: info['hydrophobic']["cluster"] for protein_id, info in feature_dict.items()}
		template_name = args.template_name  # User define P33557
		## User define
		template_active_architecture = args.template_active_res.split(",")
		# template_active_architecture = "33,35,37,64,66,91,93,97,99,106,108,115,116,118,142,146,147,148,154,156,158,191,197,199,200"

		## calculate alignment USING align_structure function
		alignment_file = align_structure(structure_folder)
		# alignment_file = r"D:\subject\active\1-qProtein\code\test\usalign_test.fasta"
		local_hydrophobic_file = os.path.join(wd, "local_hydrophobic_feature.csv")
		local_aa_file = os.path.join(wd, "local_aa_feature.csv")
		local.run(template_name, alignment_file,
		          template_active_architecture, hydrophobic_feature,
		          structure_folder, local_hydrophobic_file, local_aa_file,
		          active_edge_dist=int(args.dist1), intermediate_edge_dist=int(args.dist2))

	t2 = time.time()
	using_time = t2-t1
	logger.info(f"qProtein analysis finished in {int(using_time)} seconds!")
	
main()
