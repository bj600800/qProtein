"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2024/11/27

# Description: hydrophobic cluster analysis for typical hydrophobic residues, val, ile, leu
# ------------------------------------------------------------------------------
"""
import warnings

import biotite.structure as struc
import networkx as nx
import numpy as np
from biotite.structure.info import vdw_radius_single

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


def resi_area():
    r_vdw = vdw_radius_single("C")
    pi = np.pi
    ile_area = 4 * 4 * pi * r_vdw ** 2
    leu_area = 4 * 4 * pi * r_vdw ** 2
    val_area = 3 * 4 * pi * r_vdw ** 2
    area = {'ILE': ile_area, 'LEU': leu_area, 'VAL': val_area}
    return area


def detect_hydrophobic_cluster(structure, bias=1.1):
    r_vdw = vdw_radius_single("C")
    hydropho_dist = r_vdw * 2 + bias  # hydropho_dist==4.5

    hydrophobic_mask = np.isin(structure.res_name, ["ILE", "LEU", "VAL"]) & np.isin(structure.atom_name,
                                                                                    ["CB", "CG1", "CG2", "CD1", "CD2"])
    cell_list = struc.CellList(
        structure,
        cell_size=hydropho_dist,
        selection=hydrophobic_mask
    )

    res_pairs = set()
    for atom_idx in np.where(hydrophobic_mask)[0]:
        target_res_id = structure[atom_idx].res_id
        target_res_name = structure[atom_idx].res_name
        atoms_in_cellist = cell_list.get_atoms(coord=structure.coord[atom_idx], radius=hydropho_dist)
        potential_bond_partner_indices = [idx for idx in atoms_in_cellist if structure[idx].res_id != target_res_id]
        for potential_atom_idx in potential_bond_partner_indices:
            potential_res_id = structure[potential_atom_idx].res_id
            potential_res_name = structure[potential_atom_idx].res_name
            res_pairs.add(tuple(sorted([(target_res_id, target_res_name),
                                        (potential_res_id, potential_res_name)])))
    return list(sorted(res_pairs))


def create_hydrophobic_graph(res_pairs):
    graph = nx.Graph()
    res_list = sorted([i for res_pair in res_pairs for i in res_pair])

    for res in res_list:
        graph.add_node(res[0], res_name=res[1])

    for res1, res2 in res_pairs:
        res1_id, res2_id = res1[0], res2[0]
        graph.add_edge(res1_id, res2_id)

    return graph


def calculate_weighted_sum_of_clusters(G):
    residue_area = resi_area()
    connected_components = nx.connected_components(G)
    cluster = {}
    for idx, component in enumerate(connected_components):
        if len(component) > 1:
            component_subgraph = G.subgraph(component)

            nodes = component_subgraph.nodes()
            res_ids = [int(i) for i in nodes]

            node_labels = nx.get_node_attributes(G, 'res_name')
            res_names = [node_labels[node] for node in nodes]
            sum_area = '{:.2f}'.format(sum([residue_area[i] for i in res_names]))
            cluster[f'cluster_{idx}'] = [sum_area, res_ids]
    return cluster


def run(atom_array):
    res_pairs = detect_hydrophobic_cluster(atom_array)
    Graph = create_hydrophobic_graph(res_pairs)
    cluster = calculate_weighted_sum_of_clusters(Graph)
    return cluster