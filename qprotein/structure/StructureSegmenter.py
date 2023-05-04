"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/03/14

# Description: query structure segmentation
# ------------------------------------------------------------------------------
"""

import tempfile
import os
from numpy import mean
from pdbecif.mmcif_tools import MMCIF2Dict

from qprotein.utilities import logger
logger = logger.setup_log(name=__name__)


class Segmenter:
    def __init__(self, cif_string, query_name, subject_start_residue, subject_end_residue, query_start_residue,
                 query_end_residue, query_match_length, query_length, write_cif_path_kw, uniprot_class):
        self.cif_dict = None
        self.temp_cif_path = None
        self.subj_structure_length = None
        self.query_name = query_name
        self.cif_string = cif_string
        self.subject_start_residue = int(subject_start_residue)
        self.subject_end_residue = int(subject_end_residue)
        self.query_start_residue = int(query_start_residue)
        self.query_end_residue = int(query_end_residue)
        self.query_match_length = int(query_match_length)
        self.query_length = int(query_length)
        self.write_cif_path_kw = write_cif_path_kw
        self.uniprot_class = uniprot_class
        self.dict_to_file()

    def dict_to_file(self):
        mmcif_dict = MMCIF2Dict()
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        temp_file.write(self.cif_string.encode())
        temp_file.flush()
        self.temp_cif_path = temp_file.name
        self.cif_dict = mmcif_dict.parse(self.temp_cif_path)
        temp_file.close()

    @staticmethod
    def get_plddt_obj(structure_obj):
        bfactor_list = []
        for model in structure_obj:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        bfactor = atom.get_bfactor()
                        if bfactor is not None:
                            bfactor_list.append(bfactor)
        mean_plddt = sum(bfactor_list) / len(bfactor_list)
        return mean_plddt

    def get_plddt_dict(self, mmcif_dict, cif_id):
        bfactor = [float(i) for i in mmcif_dict[cif_id]['_atom_site']['B_iso_or_equiv']]
        avg_plddt = mean(bfactor)
        return avg_plddt

    def left_unmatched(self):
        # print("###left_unmatched")
        add_residue = self.query_start_residue - 1
        if self.subject_start_residue - add_residue > 0:
            start_residue_number = self.subject_start_residue - add_residue
        else:
            start_residue_number = 1
        return start_residue_number

    def right_unmatched(self):
        # print("###right_unmatched")
        add_residue = self.query_length - self.query_end_residue
        if self.subject_end_residue + add_residue <= self.subj_structure_length:
            end_residue_number = self.subject_end_residue + add_residue
        else:
            end_residue_number = self.subj_structure_length
        return end_residue_number

    def get_start_end_residue(self):
        start_residue_number = self.subject_start_residue
        end_residue_number = self.subject_end_residue
        # start-part missing
        if self.query_start_residue > 1:
            start_residue_number = self.left_unmatched()

        # right-part missing
        if self.query_end_residue < self.query_length:
            end_residue_number = self.right_unmatched()
        return start_residue_number, end_residue_number

    @staticmethod
    def get_cif_dict(cif_path):
        mmcif_dict = MMCIF2Dict()
        cif_dict = mmcif_dict.parse(cif_path)
        return cif_dict

    def get_atom_site(self, start_residue_number, end_residue_number):
        cif_dict = self.get_cif_dict(self.temp_cif_path)
        start_num = []
        end_num = []
        cif_id = list(cif_dict.keys())[0]
        for num, res_id in enumerate(cif_dict[cif_id]['_atom_site']['label_seq_id']):
            if int(res_id) == start_residue_number:
                start_num.append(num)
            if int(res_id) == end_residue_number:
                end_num.append(num)
        return min(start_num), max(end_num)

    def get_cif_info(self):
        cif_dict = self.get_cif_dict(self.temp_cif_path)
        cif_id = list(cif_dict.keys())[0]
        self.subj_structure_length = int(cif_dict[cif_id]['_atom_site']['label_seq_id'][-1])
        entry = cif_dict[cif_id]['_entry']
        return cif_dict, cif_id, entry

    @staticmethod
    def extract_structure(cif_dict, cif_id, entry, start_residue_number, start_atom, end_atom):
        atom_site_dict = {k: v for k, v in cif_dict[cif_id]['_atom_site'].items()}
        output_cif = {cif_id: {'_entry': entry, '_atom_site': atom_site_dict}}
        # extract structure
        for k, v in output_cif[cif_id]['_atom_site'].items():
            output_cif[cif_id]['_atom_site'][k] = v[start_atom:end_atom + 1]
        for k, v in output_cif[cif_id]['_atom_site'].items():
            if k == 'id':
                output_cif[cif_id]['_atom_site'][k] = [str(i) for i in range(1, len(v)+1)]
            if k in ('label_seq_id', 'auth_seq_id'):
                output_cif[cif_id]['_atom_site'][k] = [int(i) - start_residue_number + 1 for i in v]
        return output_cif

    def get_output_structure(self):
        cif_dict, cif_id, entry = self.get_cif_info()
        start_residue_number, end_residue_number = self.get_start_end_residue()
        start_atom, end_atom = self.get_atom_site(start_residue_number=start_residue_number, end_residue_number=end_residue_number)
        output_cif = self.extract_structure(cif_dict=cif_dict, cif_id=cif_id, entry=entry, start_residue_number=start_residue_number,
                                            start_atom=start_atom, end_atom=end_atom)
        return output_cif, cif_id

    def run(self):
        output_cif, cif_id = self.get_output_structure()
        float_plddt = float(self.get_plddt_dict(output_cif, cif_id))
        mean_plddt = round(float_plddt, 1)
        os.remove(self.temp_cif_path)
        return output_cif, mean_plddt
