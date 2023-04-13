#!/usr/bin/python3

from __future__ import print_function
import os.path
import pandas as pd
import mrcfile

from proteinpy.structure import Structure
from proteinpy.contact import Contact
from proteinpy.network import Network
from proteinpy.electron_density import ElectronDensity

class Parser:
    """Parse a PDB/DSSP/CCP4 file and returns a DataFrame object."""
    def __init__(self, pdb_id, pdb_only=False):
        """Create a Parser object.
        The Parser call a number of functions.
        """
        self._pdb_id = pdb_id
        self._pdb_df = pd.DataFrame(columns=['PDB','CHAIN_ID','RES_SEQ_N','RES_NAME',\
                                             'ATOM_NAME','COOR_X','COOR_Y','COOR_Z','DSSP'])
        self.electrondensity = ElectronDensity()
        
        if pdb_only:
            try:
                self.fetch_pdb()
            except:
                pass
        else:
            try:
                self.fetch_pdb        
                self.fetch_dssp()
                self.fetch_ccp4()
            except:
                pass

        self.structure = Structure(self._pdb_df)
        self.contact = Contact(self.structure)
        self.network = Network(self.contact)
        # self.electrondensity = ElectronDensity()


    def fetch_pdb(self):
        """Parse a PDB file and returns a DataFrame object."""
        pdb_file = '../PDB/{}.pdb'.format(self._pdb_id)

        global res_names, atom_names
        # A list of dataframe column names.
        column_names = ["PDB", "ATOM", "ATOM_SN", "ATOM_NAME", "ALT_LOC_ID", \
                        "RES_NAME", "CHAIN_ID", "RES_SEQ_N", "CHAIN_RES_N", \
                        "INSERT_RES", "COOR_X", "COOR_Y", "COOR_Z", "OCCUP", \
                        "B_FACTOR", "SEGMENT_ID", "ELEMENT_SYMBOL"]

        retcol_names = ['PDB', 'CHAIN_ID', 'RES_SEQ_N', 'RES_NAME', 'ATOM_NAME', 'COOR_X',' COOR_Y', 'COOR_Z']

        res_names = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', \
                     'GLN', 'GLU', 'GLY', 'HIS', 'ILE', \
                     'LEU', 'LYS', 'MET', 'PHE', 'PRO', \
                     'SER', 'THR', 'TRP', 'TYR', 'VAL', \
                     'HOH']
        atom_names = ['N', 'C', 'CA', 'O', 'CB', 'OXT', \
                      'SD', 'CZ', 'CH2', 'CG1', 'NZ', 'ND2', \
                      'OD1', 'CD1', 'CE1', 'CE', 'CD2', 'CZ3', \
                      'SG', 'OE2', 'NE1', 'NE', 'CE2', 'CG', \
                      'NH2', 'OE1', 'CD', 'ND1', 'OH', 'NE2', \
                      'NH1', 'CG2', 'OD2', 'CZ2', 'CE3', 'OG', 'OG1']

        if not os.path.exists(pdb_file):
            print("ERROR! CANNOT FIND PDB FILE.")
            return self

        l = []
        with open(pdb_file, 'r') as f:
            for line in f.readlines():
                if (line[:6].strip() in ["ATOM","HETATM"]) and (line[17:20].strip() in res_names) and (line[12:16].strip() in atom_names): # select all atoms (including water&residue) for pdb file; only water for symmetry pdb file
                    ATOM_CODE = line[:6].strip() # e.g. "HETATM"
                    ATOM_SN = line[6:11].strip() # Atom serial number
                    ATOM_NAME = line[12:16].strip() # e.g. "O"
                    ALT_LOC_ID = line[16].strip() # Alternative location indicator
                    RES_NAME = line[17:20].strip() # e.g. "HOH"
                    CHAIN_ID = line[21].strip() # e.g. "A"
                    RES_SEQ_N = line[22:26].strip() # Residue sequence number, e.g. "499"
                    CHAIN_RES_N = CHAIN_ID+'.'+RES_SEQ_N
                    INSERT_RES = line[26]
                    COOR_X, COOR_Y, COOR_Z = line[30:38], line[38:46], line[46:54] # "-11.662, 44.101, -1.975"
                    OCCUP = line[54:60]
                    B_FACTOR = line[60:66]
                    SEGMENT_ID = line[72:76]
                    ELEMENT_SYMBOL = line[76:78]
                    
                    l.append([self._pdb_id, ATOM_CODE, ATOM_SN,ATOM_NAME,ALT_LOC_ID,RES_NAME,CHAIN_ID,RES_SEQ_N,CHAIN_RES_N,INSERT_RES,COOR_X,COOR_Y,COOR_Z,OCCUP,B_FACTOR,SEGMENT_ID,ELEMENT_SYMBOL])

        pdb_df = pd.DataFrame(l, columns=column_names)
        pdb_df['RES_SEQ_N'] = pdb_df['RES_SEQ_N'].astype('int64')
        pdb_df[["COOR_X", "COOR_Y", "COOR_Z"]] = pdb_df[["COOR_X", "COOR_Y", "COOR_Z"]].astype('float')
        pdb_df.drop_duplicates(subset=['COOR_X', 'COOR_Y', 'COOR_Z'], keep='first', inplace=True)

        df = pdb_df[['PDB','CHAIN_ID','RES_SEQ_N','RES_NAME','ATOM_NAME','COOR_X','COOR_Y','COOR_Z']]
        self._pdb_df = df
        return self
    
    def fetch_dssp(self):#, ccp4_file='../CCP4/{}.ccp4'.format(self.pdb_id)):
        """Parse a DSSP file and returns a DataFrame object."""
        dssp_file = '../DSSP/{}.dssp'.format(self._pdb_id)

        column_names = ["PDB", "RES_NAME", "CHAIN_ID", "RES_SEQ_N", "DSSP", "S_STRUCTURE"]
        res_dict = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', \
                    'Q':'GLN', 'E':'GLU', 'G':'GLY', 'H':'HIS', 'I':'ILE', \
                    'L':'LEU', 'K':'LYS', 'M':'MET', 'F':'PHE', 'P':'PRO', \
                    'S':'SER', 'T':'THR', 'W':'TRP', 'Y':'TYR', 'V':'VAL'}
        dssp_dict = {'H':'Alpha helix', 'B':'Beta bridge', \
                     'E':'Strand', 'G':'Helix-3', \
                     'I':'Helix-5', 'T':'Turn', \
                     'S':'Bend', ' ':'Loop', \
                     'NIL':'NIL'}
        ss_dict = {column: [] for column in column_names}

        if not os.path.exists(dssp_file):
            print("ERROR! CANNOT FIND DSSP FILE.")
            return self

        f = open(dssp_file, 'r')
        datastarted = False
        for line in f.readlines():
            if line.split()[0]=="#" and line.split()[1]=="RESIDUE":
                datastarted = True
            if datastarted and line.split()[0]!="#" and (line[13:14] in res_dict.keys() or line[13:14].islower()):
                RES_NAME = res_dict[line[13:14]] if line[13:14].isupper() else 'CYS'
                CHAIN_ID = line[11:12]
                RES_SEQ_N = (line[6:10])
                RESI_N = CHAIN_ID + "."+ (RES_SEQ_N)
                S_STRUCTURE = dssp_dict[line[16:17]]
                ss_dict['PDB'].append(self._pdb_id)
                ss_dict["RES_NAME"].append(RES_NAME)
                ss_dict["CHAIN_ID"].append(CHAIN_ID)
                ss_dict["RES_SEQ_N"].append(RES_SEQ_N)
                ss_dict["DSSP"].append(line[16:17])
                ss_dict["S_STRUCTURE"].append(S_STRUCTURE)
        ss_df = pd.DataFrame(ss_dict, columns=column_names)
        ss_df['RES_SEQ_N'] = ss_df['RES_SEQ_N'].astype('int64')

        self._pdb_df = self._pdb_df.merge(ss_df.iloc[:, :-1], how='left', on=['PDB','RES_NAME','CHAIN_ID','RES_SEQ_N',]).fillna('NA')
        return self


    def fetch_ccp4(self):#, ccp4_file='../CCP4/{}.ccp4'.format(self.pdb_id)):
        self.electrondensity = self.electrondensity.fetch_data(pdb_id=self._pdb_id)
        return self


