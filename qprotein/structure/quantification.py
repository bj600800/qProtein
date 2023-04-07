#!/usr/bin/env python           
# -*- coding:utf-8 -*-          
# @Filename:    pdb_parser.py      
# @Author:      Eric Dou        
# @Time:        2021/5/27 9:02 

"""Function: Information extraction from .pdb file, i.e. sequence, dss, etc."""
import os
import csv
from Bio.PDB import MMCIFParser
from Bio.PDB.DSSP import DSSP
import pdb2pqr
from tqdm import tqdm
import freesasa
import io

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import PDBParser, MMCIFIO


def get_dssp_cif(struct_path, dssp_path):
    parser = MMCIFParser()
    struct_file_name = os.path.split(struct_path)[1].split('.')[0]
    model = parser.get_structure(struct_file_name, struct_path)[0]
    dssp = DSSP(model, struct_path, dssp_path)
    dssp_dat = [residue_dat for residue_dat in dssp]
    return dssp_dat


def get_dssp_pdb(struct_path, dssp_path):
    p = PDBParser()
    struct_file_name = os.path.split(struct_path)[1].split('.')[0]
    structure = p.get_structure(struct_file_name, struct_path)
    model = structure[0]
    dssp = DSSP(model, struct_path, dssp_path)
    dssp_dat = [residue_dat for residue_dat in dssp]
    return dssp_dat


def get_second_struct(dssp_dat):
    second_struct = [residue_dat[2].replace('-', 'C') for residue_dat in dssp_dat]
    ss = ''.join(second_struct).replace('G', 'H').replace('I', 'H')
    return ss


def get_ss_content(ss):
    ss_content = {}
    total_length = len(ss)
    ss_content['helix'] = ss.count('H')/total_length
    ss_content['sheet'] = ss.count('E')/total_length
    ss_content['loop'] = ss.count('C') / total_length
    ss_content['turn'] = ss.count('T')/total_length
    ss_content['bend'] = ss.count('S')/total_length
    ss_content['bridge'] = ss.count('B')/total_length
    return ss_content


def get_backbone_geometry(ss):
    ss_sequence = []
    i = 0
    while i < len(ss):
        c = ss[i]
        count = 1
        while i + count < len(ss) and ss[i + count] == c:
            count += 1
        ss_sequence.append(c + c * (count - 1))
        i += count


def get_Hbond(dssp_dat):
    """
    surface: Rsass > 0.05
    """
    aa_num = len(dssp_dat)
    hbond_num = 0
    hbond_energy = 0
    hbond_cols = [6, 8, 10, 12]
    hbond_energy_cols = [7, 9, 11, 13]
    for residue_dat in dssp_dat:
        for i, col in enumerate(hbond_cols):
            if residue_dat[col] != 0:
                hbond_num += 1
                hbond_energy += residue_dat[hbond_energy_cols[i]]
    hbond = {'hbond_density': hbond_num/aa_num, 'hbond_avg_energy': hbond_energy/aa_num}
    return hbond


def analyze_surface(dssp_dat):
    """
    surface: Rsass > 0.05
    the number of polar residues, non-polar residues, positive charge residues, negative charge residues.
    """
    # Define residue groups
    residue_groups = {
        'polar': ('S', 'T', 'Y', 'N', 'Q', 'C', 'K', 'R', 'H', 'D', 'E'),
        'apolar': ('G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P'),
        'positive': ('K', 'R', 'H'),
        'negative': ('D', 'E')
        }
    # Filter surface residues
    surface_residues = [(dat[0], dat[1]) for dat in dssp_dat if dat[3] > 0.05]
    total_surface_residue = len(surface_residues)

    # Calculate proportions of different residue groups
    surface = {}
    for group, residues in residue_groups.items():
        count_group = len([res_site for res_site, res_type in surface_residues if res_type in residues])
        surface[group] = count_group / total_surface_residue
    return surface


def cif2pdb(struct_path):
    struct_name = os.path.splitext(os.path.basename(struct_path))[0]
    parser = MMCIFParser()
    structure = parser.get_structure(struct_name, struct_path)
    pdb_output = io.StringIO()
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(pdb_output)
    pdb_file = io.StringIO(pdb_output.getvalue())
    return pdb_file


def get_sasa(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("pdb", pdb_file)
    result, sasa = freesasa.calcBioPDB(structure)
    # sasa-> {'Polar': 4657.344621590802, 'Apolar': 5662.672941379591}
    return sasa


def get_electrostatics():
    pass


def process_structure(struct_path):
    """
    Iterate every structure for process.
    """
    dssp = r'D:\subject\active\1-qProtein\code\qprotein\structure\tools\dssp\dssp-3.0.0-win32.exe'
    if struct_path.split('.')[-1] == 'cif':
        dssp_dat = get_dssp_cif(struct_path, dssp)
        pdb_file = cif2pdb(struct_path)
        sasa = get_sasa(pdb_file)
    else:
        dssp_dat = get_dssp_pdb(struct_path, dssp)
        sasa = get_sasa(struct_path)
    ss_content = get_ss_content(get_second_struct(dssp_dat))
    hbond = get_Hbond(dssp_dat)
    surface = analyze_surface(dssp_dat)
    return {'ss_content': ss_content, 'hbond': hbond, 'surface': surface, 'sasa': sasa}


def write_to_csv(struct_dir, results_output):
    struct_path_list = [os.path.join(struct_dir, i) for i in os.listdir(struct_dir)]
    total = len(struct_path_list)
    with open(results_output, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['', 'helix', 'sheet', 'loop', 'turn', 'bend', 'bridge', 'hbond_density',
                         'hbond_avg_energy', 'apolar', 'polar', 'positive', 'negative',
                         'polar_area', 'apolar_area'])
    for struct_path in tqdm(struct_path_list, desc='Progress', unit='step'):
        struct_name = os.path.split(struct_path)[1].split('.')[0]
        results = process_structure(struct_path)
        writer.writerow([struct_name]+[str(value) for subdict in results.values() for value in subdict.values()])


if __name__ == '__main__':
    struct_dir = r'D:\subject\active\1-qProtein\data\tibet\structure90'
    results_output = r'D:\subject\active\1-qProtein\data\tibet\structure90_results.csv'
    write_to_csv(struct_dir, results_output)











