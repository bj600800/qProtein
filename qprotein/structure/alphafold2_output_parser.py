#!/usr/bin/env python           
# -*- coding:utf-8 -*-          
# @Filename:    service.py      
# @Author:      Eric Dou        
# @Time:        2022/4/18 8:52 

"""for service, batch analysis for alphafold"""

import argparse
import json
import os
import shutil
import warnings

import pandas as pd
from Bio import BiopythonWarning, SeqIO
from numpy import mean

warnings.simplefilter('ignore', BiopythonWarning)


class Extractor:
    def __str__(self):
        # TODO: print basic info through stdout.
        structures_num = len(self.struc_dirs)
        average_plddt = mean([int(i) for i in self.extract_plddt().values()])
        str_info = f'Found {structures_num} predicted protein structures ' \
                   f'with an average pLDDT of {average_plddt}'
        return str_info

    def __init__(self, output_dir):

        # alphafold output dir
        self.info = {}
        self.output_dir = output_dir

        # list for every structure dir absPath
        self.struc_dirs = [os.path.join(self.output_dir, i) for i in os.listdir(self.output_dir)
                           if os.path.exists(os.path.join(self.output_dir, i, 'ranking_debug.json'))]
        # Create dictionary to store results
        self.make_dir()

    def make_dir(self):
        if not os.path.exists(os.path.join(self.output_dir, 'results')):
            os.mkdir(os.path.join(self.output_dir, 'results'))
        if not os.path.exists(os.path.join(self.output_dir, 'results', 'structures')):
            os.mkdir(os.path.join(self.output_dir, 'results', 'structures'))

    def copy_PDB(self):
        for path in self.struc_dirs:
            pdb_path = os.path.join(path, 'ranked_0.pdb')
            src_file = open(pdb_path, 'rb')
            tar_file = open(os.path.join(self.output_dir, 'results', 'structures',
                                         os.path.basename(path) + '.pdb'
                                         ), 'wb')
            shutil.copyfileobj(src_file, tar_file)

    def get_len(self):
        length = {}
        for path in self.struc_dirs:
            pdb_path = os.path.join(path, 'ranked_0.pdb')
            for record in SeqIO.parse(pdb_path, "pdb-atom"):
                length[os.path.basename(path)] = len(record.seq)
        return length

    def extract_plddt(self):
        # key = ncbi id, value = ranking score
        plddt = {}
        for pdb_file in self.struc_dirs:
            with open(os.path.join(pdb_file, 'ranking_debug.json'), 'r') as jsonfile:
                score = json.load(jsonfile)
                plddt[os.path.basename(pdb_file)] = max([v for k, v in score['plddts'].items()])
        return plddt

    def sum_all(self):
        self.info = {
            'Sequence length': self.get_len(),
            'pLDDT from local': self.extract_plddt(),
            }

        df = pd.DataFrame.from_dict(self.info)
        csv_path = os.path.join(self.output_dir, 'results', 'data.csv')
        df.to_csv(csv_path, sep=',')

        return print(f'Data has been successfully written to {csv_path}')

    def run(self):
        self.copy_PDB()
        self.sum_all()


parser = argparse.ArgumentParser(description='Extractor for AlphaFold2 predictions')

parser.add_argument('--output', help='output directory for AlphaFold2 prediction', required=True)

args = parser.parse_args()

output = args.output
ala = Extractor(output)
print(ala)
ala.run()
