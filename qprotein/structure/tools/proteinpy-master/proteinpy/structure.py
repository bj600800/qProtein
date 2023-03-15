#!/usr/bin/python3

from __future__ import print_function
import os.path
import pandas as pd

class Structure:
	"""Parse a PDB/DSSP/CCP4 file and returns a DataFrame object."""
	def __init__(self, _pdb_df):
		self.allatoms = _pdb_df
		self.resatoms = self.get_resatoms()
		self.waters = self.get_waters()
		self.buriedwaters = self.get_buriedwaters()

	def get_resatoms(self):

		pdb_df = self.allatoms

		global backbone_atoms, sidechain_atoms
		backbone_atoms = ['N', 'C', 'CA', 'O', 'CB', 'OXT']
		sidechain_atoms = ['SD', 'CZ', 'CH2', 'CG1', 'NZ', \
						   'ND2', 'OD1', 'CD1', 'CE1', 'CE', \
						   'CD2', 'CZ3', 'SG', 'OE2', 'NE1', \
						   'NE', 'CE2', 'CG', 'NH2', 'OE1', \
						   'CD', 'ND1', 'OH', 'NE2', 'NH1', \
						   'CG2', 'OD2', 'CZ2', 'CE3', 'OG', 'OG1']
		res_names = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', \
					 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', \
					 'LEU', 'LYS', 'MET', 'PHE', 'PRO', \
					 'SER', 'THR', 'TRP', 'TYR', 'VAL']

		resatoms = pdb_df.loc[(pdb_df["ATOM_NAME"].isin(backbone_atoms+sidechain_atoms))&\
							  (pdb_df["RES_NAME"].isin(res_names))]\
						 .drop_duplicates().copy().reset_index(drop=True)

		return resatoms

	def get_waters(self):

		pdb_df = self.allatoms
		waters = pdb_df.loc[(pdb_df["ATOM_NAME"]=='O')&\
							(pdb_df["RES_NAME"]=='HOH')]\
					   .drop_duplicates().copy().reset_index(drop=True)
		return waters

	def get_buriedwaters(self):
		pass

	def get_coordinate(self, atom):
		_pdb_id, chain_id, res_seq_n, res_name, atom_name = atom.split('.')
		pdb_df = self.allatoms

		return pdb_df.loc[(pdb_df.CHAIN_ID==chain_id)&\
						  (pdb_df.RES_SEQ_N==int(res_seq_n))&\
						  (pdb_df.RES_NAME==res_name)&\
						  (pdb_df.ATOM_NAME==atom_name), ['COOR_X','COOR_Y','COOR_Z']].values[0]
