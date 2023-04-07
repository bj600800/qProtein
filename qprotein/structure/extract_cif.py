"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/14

# Description: 
# ------------------------------------------------------------------------------
"""
import os
from numpy import mean

from pdbecif.mmcif_tools import MMCIF2Dict
from pdbecif.mmcif_io import CifFileWriter


def read_cif(cif_path):
    mmcif_dict = MMCIF2Dict()
    cif_dict = mmcif_dict.parse(cif_path)
    return cif_dict

"""
group_PDB 5021
id 5021
type_symbol 5021
label_atom_id 5021
label_alt_id 5021
label_comp_id 5021
label_asym_id 5021
label_entity_id 5021
label_seq_id 5021 -> residue number
pdbx_PDB_ins_code 5021
Cartn_x 5021
Cartn_y 5021
Cartn_z 5021
occupancy 5021
B_iso_or_equiv 5021
pdbx_formal_charge 5021
auth_seq_id 5021
auth_comp_id 5021
auth_asym_id 5021
auth_atom_id 5021
pdbx_PDB_model_num 5021
pdbx_sifts_xref_db_acc 5021
pdbx_sifts_xref_db_name 5021
pdbx_sifts_xref_db_num 5021 -> residue number
pdbx_sifts_xref_db_res 5021
"""


def get_position_num(cif_dict, start_res, end_res):
    start_num = []
    end_num = []
    cif_id = list(cif_dict.keys())[0]
    for num, res_id in enumerate(cif_dict[cif_id]['_atom_site']['label_seq_id']):
        if res_id == start_res:
            start_num.append(num)
        if res_id == end_res:
            end_num.append(num)

    return (min(start_num), max(end_num))


def edit_cif(start_res, atom_num, cif_dict):
    """
    _entity_poly_seq, _ma_qa_metric_local
    """
    cif_id = list(cif_dict.keys())[0]
    mean_bfactor = mean([float(i) for i in cif_dict[cif_id]['_atom_site']['B_iso_or_equiv']])
    if mean_bfactor < 90:
        return
    entry = cif_dict[cif_id]['_entry']
    start_atom = atom_num[0]
    end_atom = atom_num[1]

    atom_site_dict = {k: v for k, v in cif_dict[cif_id]['_atom_site'].items()}
    output_cif = {cif_id: {'_entry': entry, '_atom_site': atom_site_dict}}

    # extract structure
    for k, v in output_cif[cif_id]['_atom_site'].items():
        output_cif[cif_id]['_atom_site'][k] = v[start_atom:end_atom+1]

    # reorder residue number
    for k, v in output_cif[cif_id]['_atom_site'].items():
        if k in ('label_seq_id', 'auth_seq_id'):
            output_cif[cif_id]['_atom_site'][k] = [int(i) - start_res + 1 for i in v]
    return output_cif


def write_cif(cif_dict, write_path):
    if cif_dict:
        writer = CifFileWriter(write_path)
        writer.write(cif_dict)


if __name__ == '__main__':
    dir_path = r'D:\subject\active\1-qProtein\data\tibet\sprot90_90'
    write_dir = r'D:\subject\active\1-qProtein\data\tibet\sprot90_90_segment'
    cath_file = r'D:\subject\active\1-qProtein\data\tibet\ident90_from_cif.out'
    if not os.path.exists(write_dir):
        os.mkdir(write_dir)
    with open(cath_file, 'r') as rf:
        content = rf.readlines()
    for cif in os.listdir(dir_path):
        for line in content:
            line = line.split('\t')
            name = line[0]
            start_res = line[6]
            end_res = line[7]
            if int(line[3]) > 280:
                cif_name = cif.split('.')[0]
                cif_path = os.path.join(dir_path, cif)
                if cif_name == name:
                    cif_dict = read_cif(cif_path)
                    atom_num = get_position_num(cif_dict, start_res, end_res)
                    cif_dict = edit_cif(int(start_res), atom_num, cif_dict)
                    write_path = os.path.join(write_dir, cif_name+'.cif')
                    print(write_path)
                    write_cif(cif_dict, write_path)
