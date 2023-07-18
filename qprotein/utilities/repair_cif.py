"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/07/11

# Description: 
# ------------------------------------------------------------------------------
"""
import os

from pdbecif.mmcif_tools import MMCIF2Dict
from pdbecif.mmcif_io import CifFileWriter
from Bio.SeqUtils import IUPACData

def get_cif_info(cif_path):
	mmcif_dict = MMCIF2Dict()
	cif_dict = mmcif_dict.parse(cif_path)
	cif_id = list(cif_dict.keys())[0]
	entry = cif_dict[cif_id]['_entry']
	return cif_dict, cif_id, entry

def repair_cif(cif_dict, cif_id, entry):
	atom_site_dict = {k: v for k, v in cif_dict[cif_id]['_atom_site'].items()}
	output_cif = {cif_id: {'_entry': entry, '_atom_site': atom_site_dict}}
	# extract structure
	output_cif[cif_id]['_atom_site']['auth_seq_id'] = output_cif[cif_id]['_atom_site']['label_seq_id']
	output_cif[cif_id]['_atom_site']['auth_comp_id'] = output_cif[cif_id]['_atom_site']['label_comp_id']
	output_cif[cif_id]['_atom_site']['auth_atom_id'] = output_cif[cif_id]['_atom_site']['label_atom_id']
	output_cif[cif_id]['_atom_site']['pdbx_sifts_xref_db_acc'] = ['protein']*len(output_cif[cif_id]['_atom_site']['label_atom_id'])
	output_cif[cif_id]['_atom_site']['pdbx_sifts_xref_db_name'] = ['UNP']*len(output_cif[cif_id]['_atom_site']['label_atom_id'])
	output_cif[cif_id]['_atom_site']['pdbx_sifts_xref_db_num'] = output_cif[cif_id]['_atom_site']['label_seq_id']
	output_cif[cif_id]['_atom_site']['pdbx_sifts_xref_db_res'] = [IUPACData.protein_letters_3to1[i[0]+i[1:].lower()] for i in output_cif[cif_id]['_atom_site']['label_comp_id']]
	return output_cif

def write_cif(write_cif_path, output_cif):
	writer = CifFileWriter(write_cif_path)
	writer.write(output_cif)

def run(src_cif, out_cif):
	cif_dict, cif_id, entry = get_cif_info(src_cif)
	output_cif = repair_cif(cif_dict, cif_id, entry)
	write_cif(out_cif, output_cif)

if __name__ == '__main__':
	src_dir = r"D:\subject\active\1-qProtein\data\enzymes\cazy_endo-1_4-beta-xylanase\classification\src_dir"
	out_dir = r"D:\subject\active\1-qProtein\data\enzymes\cazy_endo-1_4-beta-xylanase\classification\out_dir"
	src_list = os.listdir(src_dir)
	for file in src_list:
		src_cif = os.path.join(src_dir, file)
		out_cif = os.path.join(out_dir, file)
		run(src_cif, out_cif)