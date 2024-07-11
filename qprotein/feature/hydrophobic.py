"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/07/25

# Description: Service of hydrophobic cluster detection.
# ------------------------------------------------------------------------------
"""

import numpy as np
import networkx as nx
import biotite.structure as struc
import biotite.structure.io as strucio
from biotite.structure.info import vdw_radius_single
import warnings

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
    # 筛选满足要求的残基和原子, 原子序号从0开始。
    hydrophobic_mask = np.isin(structure.res_name, ["ILE", "LEU", "VAL"]) & np.isin(structure.atom_name,
                                                                                    ["CB", "CG1", "CG2", "CD1", "CD2"])
    cell_list = struc.CellList(
        structure,
        cell_size=hydropho_dist,
        selection=hydrophobic_mask
    )
    # 残基对之间的相互作用原子对
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
    # 使用图来描述蛋白质疏水作用，每一个连通图是一个簇
    graph = nx.Graph()
    res_list = sorted([i for res_pair in res_pairs for i in res_pair])
    # 添加所有残基节点
    for res in res_list:
        graph.add_node(res[0], res_name=res[1])
    # 添加相互作用边
    for res1, res2 in res_pairs:
        res1_id, res2_id = res1[0], res2[0]
        graph.add_edge(res1_id, res2_id)
    return graph


def calculate_weighted_sum_of_clusters(G):
    residue_area = resi_area()
    connected_components = nx.connected_components(G)  # 获取所有连通分量
    weighted_sums = []  # 存储每个连通图的权重之和
    cluster = []
    for component in connected_components:
        if len(component) > 1:
            component_subgraph = G.subgraph(component)
            nodes = component_subgraph.nodes()
            
            node_labels = nx.get_node_attributes(G, 'res_name')
            res_names = [node_labels[node] for node in nodes ]

            # print(list(nodes))
            area = sum([residue_area[res] for res in res_names])
            keys = list(nodes)
            values = [residue_area[i] for i in res_names]
            pairs = zip(keys, values)
            cluster.append({key: value for key, value in pairs})
            # print(area)
            weighted_sums.append(area)
    # for res in cluster:
    #     res_str = '+'.join([str(i) for i in res])
    #     print(res_str)
    # input()
    weighted_sum = sum(weighted_sums)
    max_cluster = max(weighted_sums)
    # sorted_unique_lst = sorted(set(weighted_sums), reverse=True)
    # print("sum:", round(weighted_sum, 2))
    return round(weighted_sum, 2), round(max_cluster, 2), cluster


def analyze(atom_array):
    hydrophobic_cluster = {}
    length = len(set(atom_array.res_id))
    res_pairs = detect_hydrophobic_cluster(atom_array)
    Graph = create_hydrophobic_graph(res_pairs)
    weighted_sum, max_cluster, cluster_res = calculate_weighted_sum_of_clusters(Graph)
    # hydrophobic_cluster['graph'] = Graph
    hydrophobic_cluster['sum_area'] = weighted_sum/length
    hydrophobic_cluster['max_area'] = max_cluster
    hydrophobic_cluster['cluster'] = cluster_res
    return hydrophobic_cluster

if __name__ == '__main__':
	struc_path = r"C:\Users\bj600\Desktop\GI24158867.pdb"
	structure = strucio.load_structure(struc_path)
	analyze(structure)