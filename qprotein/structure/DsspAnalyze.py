"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/07/12

# Description: 
# ------------------------------------------------------------------------------
"""
import os
import io
import csv

from tqdm import tqdm
import threading
import multiprocessing
from queue import Queue

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.DSSP import DSSP
from Bio import BiopythonWarning
import warnings
warnings.simplefilter("ignore", BiopythonWarning)


def get_dssp_cif(struct_path):
    parser = MMCIFParser()
    struct_file_name = os.path.split(struct_path)[1].split('.')[0]
    model = parser.get_structure(struct_file_name, struct_path)[0]
    # dssp = DSSP(model, struct_path, dssp="mkdssp")
    dssp = DSSP(model, struct_path, dssp=r"D:\subject\active\1-qProtein\tools\dssp\dssp-3.0.0-win32.exe")
    dssp_dat = [residue_dat for residue_dat in dssp]
    return dssp_dat


def get_dssp_pdb(struct_path):
    p = PDBParser()
    struct_file_name = os.path.split(struct_path)[1].split('.')[0]
    structure = p.get_structure(struct_file_name, struct_path)
    model = structure[0]
    # dssp = DSSP(model, struct_path, dssp="mkdssp")
    dssp = DSSP(model, struct_path, dssp=r"D:\subject\active\1-qProtein\tools\dssp\dssp-3.0.0-win32.exe")
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


def run(struct_path):
    """
    Iterate every structure for process.
    """
    dssp_dat = ''
    sasa = ''
    if struct_path.split('.')[-1] == 'cif':
        dssp_dat = get_dssp_cif(struct_path)
        pdb_file = cif2pdb(struct_path)
        sasa = get_sasa(pdb_file)
    elif struct_path.split('.')[-1] == 'pdb':
        dssp_dat = get_dssp_pdb(struct_path)
        sasa = get_sasa(struct_path)
    ss_content = get_ss_content(get_second_struct(dssp_dat))
    hbond = get_Hbond(dssp_dat)
    surface = analyze_surface(dssp_dat)
    return {'ss_content': ss_content, 'hbond': hbond, 'surface': surface, 'sasa': sasa}