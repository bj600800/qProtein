#!/usr/bin/python3

from __future__ import print_function
import os.path
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from proteinpy.geometry_funcs import *

class Network:
    """Parse a PDB/DSSP/CCP4 file and returns a DataFrame object."""
    def __init__(self, contact):
        """
        
        """
        self._contact = contact
        self.allnet = self.get_allnet()
        self.waternet = self.get_waternet()
        self.watertetrahedrons = self.find_watertetrahedrons()
        self.waterpentagons = self.find_watercycles(n=5)
        self.waterhexagons = self.find_watercycles(n=6)

    @classmethod
    def get_net(cls, contact_df, structure_df):
        contact_df['ATOM_1'] = contact_df['PDB_1'] + '.' + \
                            + contact_df['CHAIN_ID_1'] + '.' + \
                            + contact_df['RES_SEQ_N_1'].astype(str) + '.' + \
                            + contact_df['RES_NAME_1'] + '.' + \
                            + contact_df['ATOM_NAME_1']
        contact_df['ATOM_2'] = contact_df['PDB_2'] + '.' + \
                            + contact_df['CHAIN_ID_2'] + '.' + \
                            + contact_df['RES_SEQ_N_2'].astype(str) + '.' + \
                            + contact_df['RES_NAME_2'] + '.' + \
                            + contact_df['ATOM_NAME_2']
        contact_df['VECTOR'] = contact_df[['COOR_X_2','COOR_Y_2','COOR_Z_2']]\
                               .apply(lambda r: tuple(r), axis=1).apply(np.array)-\
                               contact_df[['COOR_X_1','COOR_Y_1','COOR_Z_1']]\
                             .apply(lambda r: tuple(r), axis=1).apply(np.array)

        structure_df['ATOM'] = structure_df['PDB'] + '.' + \
                             + structure_df['CHAIN_ID'] + '.' + \
                             + structure_df['RES_SEQ_N'].astype(str) + '.' + \
                             + structure_df['RES_NAME'] + '.' + \
                             + structure_df['ATOM_NAME']
        structure_df['COOR'] = structure_df[['COOR_X','COOR_Y','COOR_Z']]\
                             .apply(lambda r: tuple(r), axis=1).apply(np.array)
        structure_dict = structure_df.set_index('ATOM')[['RES_NAME','ATOM_NAME','COOR']].to_dict('index')

        net = nx.from_pandas_edgelist(contact_df,\
                                      source='ATOM_1',\
                                      target='ATOM_2',\
                                      edge_attr=['DIST','VECTOR'],\
                                      create_using=None)
        nx.set_node_attributes(net, structure_dict)
        return net

    def get_allnet(self):
        contact_df = self._contact.all2all.copy()
        structure_df = self._contact._structure.allatoms.copy()
        return Network.get_net(contact_df, structure_df)

    def get_waternet(self):
        contact_df = self._contact.water2all.copy()
        structure_df = self._contact._structure.allatoms.copy()
        return Network.get_net(contact_df, structure_df)

    def find_watertetrahedrons(self):

        watercenters = [node for node in self.waternet.nodes(data=True) if self.waternet.degree(node[0])==4 and node[1]['RES_NAME']=='HOH']
        waterdict = {node[0]:node[1] for node in self.waternet.nodes(data=True)}

        watertetrahedrons = []
        for center in watercenters:
            edges = self.waternet.edges(center[0], data=True)
            watertetrahedron = Tetrahedron(edges=edges, nodes=waterdict)
            # nx.set_node_attributes(watertetrahedron, waterdict)
            watertetrahedrons.append(watertetrahedron)

        return watertetrahedrons

    def find_watercycles(self, n):
        
        cycles_ = [self.waternet.subgraph(cycle) for cycle in nx.cycle_basis(self.waternet) if len(cycle)==n]
        waterdict = {node[0]:node[1] for node in self.waternet.nodes(data=True)}

        watercycles = []
        for cycle_ in cycles_:
            edges = [(node1, node2, self.waternet.get_edge_data(node1, node2)) for node1, node2 in nx.find_cycle(cycle_)]
            watercycle = Cycle(edges=edges, nodes=waterdict)
            # nx.set_node_attributes(watercycle, waterdict)
            watercycles.append(watercycle)
        
        return watercycles

    def get_vectors(self, atom1, *atoms):
        contact = self._contact
        vectors = []
        for atom in atoms:
            vectors.append(contact.get_vector(atom1, atom))
        return np.array(vectors)

class Tetrahedron(nx.DiGraph):
    def __init__(self, edges=None, nodes=None, **attr):
        super(Tetrahedron, self).__init__(edges, **attr)
        nx.set_node_attributes(self, nodes)
        self.coordinates = self.get_coordinates()
        self.vectors = self.get_vectors()
        self.angles = self.get_angles()
        self.dihedralangles = self.get_dihedralangles()
    def get_coordinates(self):
        return [node[1] for node in self.nodes.data('COOR')]
    def get_vectors(self):
        return [edge[2] for edge in self.edges.data('VECTOR')]#list(nx.get_edge_attributes(self, 'VECTOR').values())
    def get_angles(self):
        return calc_tetrahedron_angles(self.vectors)
    def get_dihedralangles(self):
        return calc_tetrahedron_dihedral_angles(self.vectors)

class Cycle(nx.DiGraph):
    def __init__(self, edges=None, nodes=None, **attr):
        super(Cycle, self).__init__(edges, **attr)
        nx.set_node_attributes(self, nodes)
        self.coordinates = self.get_coordinates()
        self.vectors = self.get_vectors()
        self.angles = self.get_angles()
        self.dihedralangles = self.get_dihedralangles()
    def get_coordinates(self):
        return [node[1] for node in self.nodes.data('COOR')]
    def get_vectors(self):
        return [edge[2] for edge in self.edges.data('VECTOR')]#list(nx.get_edge_attributes(self, 'VECTOR').values())
    def get_angles(self):
        return calc_cycle_angles(self.vectors)
    def get_dihedralangles(self):
        return calc_cycle_dihedral_angles(self.vectors)
