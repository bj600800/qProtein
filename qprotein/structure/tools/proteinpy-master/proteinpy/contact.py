#!/usr/bin/python3

from __future__ import print_function
import os.path
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree

class Contact:
    """Parse a PDB/DSSP/CCP4 file and returns a DataFrame object."""
    def __init__(self, structure):
        """
        
        """
        self._structure = structure
        self.all2all = self.get_all2all()
        self.water2water = self.get_water2water()
        self.water2res = self.get_water2res()
        self.water2all = self.get_water2all()
        self.res2res = self.get_res2res()

        # try:
        #     self.get_water2water()
        #     self.get_water2res()
        #     self.get_res2res()
        # except:
        #     pass

    def get_all2all(self):
        df = self._structure.allatoms.copy()
        coords = np.array(df[['COOR_X','COOR_Y','COOR_Z']].values.astype(np.float))
        tree = cKDTree(coords)
        pairs = list(tree.query_pairs(r=4.9))

        idx1 = list(list(zip(*pairs))[0])
        idx2 = list(list(zip(*pairs))[1])

        df1 = df.iloc[idx1, :].copy().reset_index(drop=True).add_suffix('_1')
        df2 = df.iloc[idx2, :].copy().reset_index(drop=True).add_suffix('_2')
        dist = np.linalg.norm(coords[idx1]-coords[idx2], axis=1)

        all2all_df = pd.concat([df1, df2], axis=1)
        all2all_df['DIST'] = dist
        return all2all_df

    def get_water2water(self):
        df = self.all2all
        water2water_df = df.loc[(df.RES_NAME_1=='HOH')&\
                                 (df.ATOM_NAME_1=='O')&\
                                 (df.RES_NAME_2=='HOH')&\
                                 (df.ATOM_NAME_2=='O')&\
                                 (df.DIST<3.2)].copy().reset_index(drop=True)
        return water2water_df


    def get_water2res(self):
        df = self.all2all
        water2res_df = df.loc[(df.RES_NAME_1!=df.RES_NAME_2)&\
                              (df.ATOM_NAME_1!=df.ATOM_NAME_2)&\
                              (((df.RES_NAME_1=='HOH')&(df.ATOM_NAME_1=='O'))|\
                                (df.RES_NAME_2=='HOH')&(df.ATOM_NAME_2=='O'))&\
                              (df.DIST<3.2)].copy().reset_index(drop=True)
        return water2res_df

    def get_water2all(self):
        water2water_df = self.water2water.copy()
        water2res_df = self.water2res.copy()
        return pd.concat([water2water_df,water2res_df],axis=0)

    def get_res2res(self):
        
        df = self.all2all

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


        res2res_df = df.loc[(df.RES_NAME_1.isin(res_names))&\
                            (df.ATOM_NAME_1.isin(backbone_atoms+sidechain_atoms))&\
                            (df.RES_NAME_2.isin(res_names))&\
                            (df.ATOM_NAME_2.isin(backbone_atoms+sidechain_atoms))]

        return res2res_df

    def get_vector(self, atom1, atom2):
        structure = self._structure
        return structure.get_coordinate(atom2)-structure.get_coordinate(atom1)