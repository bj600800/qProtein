"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/14

# Description: query structure segmentation
  Input: 1. full length subject structure collected from alphafold datbase, 2. subject start and end site,
         3. query start and end site, 4. query matched length, and 5. query length

  Output: segmented structure and file name.
# ------------------------------------------------------------------------------
"""

import tempfile
import os
import io
from numpy import mean
from pdbecif.mmcif_tools import MMCIF2Dict
# from Bio.PDB import Structure
from Bio.PDB.Model import Model
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.Chain import Chain
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Structure import Structure
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


class Segmenter:
    def __init__(self, cif_string, query_name, subject_start_residue, subject_end_residue, query_start_residue,
                 query_end_residue, query_match_length, query_length, write_cif_path_kw):
        self.cif_dict = None
        self.temp_cif_path = None
        self.subj_structure_length = None
        self.output_plddt = None
        self.query_name = query_name
        self.cif_string = cif_string
        self.subject_start_residue = subject_start_residue
        self.subject_end_residue = subject_end_residue
        self.query_start_residue = query_start_residue
        self.query_end_residue = query_end_residue
        self.query_match_length = query_match_length
        self.query_length = query_length
        self.write_cif_path_kw = write_cif_path_kw

        self.get_cif_dict()

    def get_cif_dict(self):
        mmcif_dict = MMCIF2Dict()
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        temp_file.write(self.cif_string.encode())
        temp_file.flush()
        self.temp_cif_path = temp_file.name
        self.cif_dict = mmcif_dict.parse(self.temp_cif_path)
        temp_file.close()

    @staticmethod
    def get_plddt(structure_obj):
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

    def full_matched(self):
        start_residue_number = self.subject_start_residue
        end_residue_number = self.subject_end_residue
        return start_residue_number, end_residue_number

    def left_unmatched(self):
        end_residue_number = self.subject_end_residue
        if self.subject_start_residue == 1:
            start_residue_number = 1
        else:
            add_residue = self.query_start_residue - 1
            start_residue_number = min(1, self.subject_start_residue - add_residue)
        return start_residue_number, end_residue_number

    def right_unmatched(self):
        start_residue_number = 1
        if self.subject_end_residue == self.subj_structure_length:
            end_residue_number = self.subj_structure_length
        else:
            add_residue = self.subject_end_residue - self.query_length
            end_residue_number = max(self.subject_end_residue, self.subject_end_residue + add_residue)
        return start_residue_number, end_residue_number

    def both_unmatched(self):
        left_add_residue = self.query_start_residue - 1
        right_add_residue = self.subject_end_residue - self.query_length

        if self.subject_end_residue == self.subj_structure_length:
            end_residue_number = self.subject_end_residue
            start_residue_number = min(1, self.subject_start_residue - left_add_residue)
        elif self.subject_start_residue == 1:
            start_residue_number = 1
            end_residue_number = max(self.subj_structure_length, self.subject_end_residue + right_add_residue)
        else:
            start_residue_number = min(1, self.subject_start_residue - left_add_residue)
            end_residue_number = max(self.subj_structure_length, self.subject_end_residue + right_add_residue)
        return start_residue_number, end_residue_number

    def get_start_end_residue(self):
        start_residue_number = 0
        end_residue_number = 0

        # query matched completely
        if self.query_match_length == self.query_length:
            start_residue_number, end_residue_number = self.full_matched()

        # left-part missing
        elif self.query_start_residue > 1 and self.query_end_residue == self.query_length:
            start_residue_number, end_residue_number = self.left_unmatched()

        # right-part missing
        elif self.query_start_residue == 1 and self.query_end_residue < self.query_length:
            start_residue_number, end_residue_number = self.right_unmatched()

        # both side missing and add pseudo structure part
        elif self.query_start_residue > 1 and self.query_end_residue < self.query_length:
            start_residue_number, end_residue_number = self.both_unmatched()
        return start_residue_number, end_residue_number

    def extract_structure(self):
        parser = MMCIFParser()
        subj_structure = parser.get_structure("subject_structure", self.temp_cif_path)
        output_subj_structure = Structure('output')
        output_subj_model = Model(0)
        output_subj_chain = Chain("A")
        for chain in subj_structure.get_chains():
            self.subj_structure_length = len([i for i in chain.get_residues()])
            position = start_residue_number, end_residue_number = self.get_start_end_residue()
            if 0 in position:
                logger.debug('BUG HERE! - '+self.write_cif_path_kw[1])
                input()
                return
            for residue in chain.get_residues():
                if start_residue_number <= residue.id[1] <= end_residue_number:
                    output_subj_chain.add(residue)
            output_subj_model.add(output_subj_chain)
        output_subj_structure.add(output_subj_model)
        mean_plddt = round(self.get_plddt(output_subj_structure), 1)
        if mean_plddt > 70:
            write_cif_path = os.path.join(self.write_cif_path_kw[0],
                                          '#'.join(self.write_cif_path_kw[1:] + ([str(mean_plddt)]))+'.cif')
            io = MMCIFIO()
            io.set_structure(output_subj_structure)
            io.save(write_cif_path)
            logger.info(f'Write structure file for {self.query_name} in {write_cif_path}.')
            return write_cif_path

    def run(self):
        write_cif_path = self.extract_structure()
        os.remove(self.temp_cif_path)
        return write_cif_path
