#!/usr/bin/env python           
# -*- coding:utf-8 -*-          
# @Filename:    pdb_parser.py      
# @Author:      Eric Dou        
# @Time:        2021/5/27 9:02 

"""Function: Information extraction from .pdb file, i.e. sequence, dss, etc."""
import os

from Bio.PDB import MMCIFParser
from Bio.PDB.DSSP import DSSP
import pdb2pqr

import freesasa
import io

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import PDBParser, MMCIFIO


def get_dssp_dat(struct_path, dssp_path):
    parser = MMCIFParser()
    struct_file_name = os.path.split(struct_path)[1].split('.')[0]
    model = parser.get_structure(struct_file_name, struct_path)[0]
    dssp = DSSP(model, struct_path, dssp_path)
    dssp_dat = [residue_dat for residue_dat in dssp]
    return dssp_dat


def get_second_struct(dssp_dat):
    second_struct = [residue_dat[2].replace('-', 'C') for residue_dat in dssp_dat]
    ss = ''.join(second_struct).replace('G', 'H').replace('I', 'H')
    return ss


def get_ss_content(ss):
    content = {}
    total_length = len(ss)
    content['helix'] = ss.count('H')/total_length
    content['sheet_content'] = ss.count('E')/total_length
    content['loop_content'] = ss.count('C') / total_length
    content['turn_content'] = ss.count('T')/total_length
    content['bend_content'] = ss.count('S')/total_length
    content['bridge_content'] = ss.count('B')/total_length
    print('ss content', content)
    return content['helix']


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
    print('aa length:', aa_num)
    hbond_num = 0
    hbond_energy = 0
    for residue_dat in dssp_dat:
        if residue_dat[6] != 0:
            hbond_num += 1
            hbond_energy += residue_dat[7]
        if residue_dat[8] != 0:
            hbond_num += 1
            hbond_energy += residue_dat[9]
        if residue_dat[10] != 0:
            hbond_num += 1
            hbond_energy += residue_dat[11]
        if residue_dat[12] != 0:
            hbond_num += 1
            hbond_energy += residue_dat[13]

    print('hbond density:', hbond_num/aa_num)


def analyze_surface(dssp_dat):
    """
    surface: Rsass > 0.05
    the number of polar residues, non-polar residues, positive charge residues, negative charge residues.
    """
    polar = ('S', 'T', 'Y', 'N', 'Q', 'C', 'K', 'R', 'H', 'D', 'E')
    positive = ('K', 'R', 'H')
    negative = ('D', 'E')
    Apolar = ('G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P')
    surface_residues = [(residue_dat[0], residue_dat[1]) for residue_dat in dssp_dat if residue_dat[3] > 0.05]
    count_surface = len(surface_residues)
    surface_dat = {
        'apolar_num': len([residue for residue in surface_residues if residue[1] in Apolar]) / count_surface,
        'polar_num': len([residue for residue in surface_residues if residue[1] in polar])/count_surface,
        'pos_num': len([residue for residue in surface_residues if residue[1] in positive])/count_surface,
        'neg_num': len([residue for residue in surface_residues if residue[1] in negative])/count_surface
        }
    print('surface_dat', surface_dat)
    return surface_dat


def cif2pdb(struct_path):
    struct_name = os.path.split(struct_path)[1].split('.')[0]
    parser = MMCIFParser()
    structure = parser.get_structure(struct_name, struct_path)
    pdb_io = PDBIO()
    pdb_output = io.StringIO()
    pdb_io.set_structure(structure)
    pdb_io.save(pdb_output)
    pdb_str = pdb_output.getvalue()
    pdb_file = io.StringIO(pdb_str)
    return pdb_file


def get_sasa(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("pdb", pdb_file)
    result, sasa_classes = freesasa.calcBioPDB(structure)
    return sasa_classes

def get_electrostatics():
    pass


def main():
    """
    Iterate every structure for process.
    """
    struct_dir = r'D:\subject\active\1-qProtein\data\tibet\ident90_segment'
    dssp = r'D:\subject\active\1-qProtein\code\qprotein\structure\tools\dssp\dssp-3.0.0-win32.exe'
    struct_path_list = [os.path.join(struct_dir, i) for i in os.listdir(struct_dir)]
    total = len(struct_path_list)
    helix = 0
    for struct_path in struct_path_list:
        print(struct_path)
        dssp_dat = get_dssp_dat(struct_path, dssp)
        ss = get_second_struct(dssp_dat)
        helix += get_ss_content(ss)
        get_Hbond(dssp_dat)
        analyze_surface(dssp_dat)
        pdb_file = cif2pdb(struct_path)
        get_sasa(pdb_file)
    print(helix/total)

main()











