"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2025/3/28

# Description: visualization for protein structure and features
# ------------------------------------------------------------------------------
"""
import os
import subprocess

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import plotly.graph_objects as go
from pyvis.network import Network
from sklearn.decomposition import PCA

from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


def draw_graph(G):
    nx.draw(G, with_labels=True, node_color='lightblue', node_size=800, edge_color='gray', font_size=12,
            font_color='black')

    plt.show()


def draw_digraph(G, protein_name):
    # create pyvis direction network
    net = Network(notebook=True, directed=True, cdn_resources='in_line')

    # add nodes and edges
    for node in G.nodes():
        color = 'blue'
        if 'sheet' in node:
            color = '#ffd449'
        elif 'helix' in node:
            color = '#d62839'
        elif 'turn' in node:
            color = '#4c956c'
        net.add_node(node, color=color, title=node, font=dict(size=20))  # 设置字号为20

    for edge in G.edges():
        net.add_edge(edge[0], edge[1])

    # layout
    net.force_atlas_2based()
    net.show(protein_name)


def draw_graph_interactive(G, output_file):
    net = Network(notebook=True, directed=False, cdn_resources='in_line')
    G = nx.relabel_nodes(G, lambda x: str(x))
    net.from_nx(G)

    max_degree = max(dict(G.degree).values()) if G.number_of_nodes() > 0 else 1

    def get_color(degree, max_degree=max_degree):
        base_color = mcolors.hex2color("#0d47a1")
        min_color = [c + (1 - c) for c in base_color]
        factor = degree / max_degree  # normalization
        color = [(1 - factor) * min_c + factor * base_c for min_c, base_c in zip(min_color, base_color)]
        return mcolors.to_hex(color)

    for node in G.nodes():
        degree = G.degree[node]
        net.get_node(node)["color"] = get_color(degree)
        net.get_node(node)["size"] = 15
        net.get_node(node)["label"] = f"residue {node} (degree: {degree})"

    net.show(output_file)
    logger.info(f"HTML graph：{output_file}")


def show_hydrocluster_pymol(pdb_file, cluster_graphs):
    """
    Show hydrophobic clusters in PyMol.
    :return:
    """
    pymol_bin = 'pymol'
    cmd = []
    for i, cluster in enumerate(cluster_graphs):
        command = f"sele Cluster {i}, res {'+'.join([str(res_id) for res_id in cluster.nodes()])}"
        print(command.replace("+", ","))
        cmd.append(command)
    pymol_cmd = '; '.join(cmd)

    # get Anaconda base environment path
    conda_base = r"/opt/anaconda3"
    pymol_path = os.path.join(conda_base, "bin")

    # add base pymol dir to the current environmental variable path
    os.environ["PATH"] = pymol_path + os.pathsep + os.environ["PATH"]

    command_args = [pymol_bin, pdb_file, '-d', pymol_cmd]
    process = subprocess.Popen(command_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    process.communicate()


def draw_ss_3D(points, group, ss_id, protein_name, fig_dir):
    """
    draw all ss in a figure.
    To confirm the correctness of the computed regression line and direction vector.

    :param points: {ss_id:[[x,y,z], [x,y,z]], ss_id:[[x,y,z], [x,y,z]]}
    :param group: {'residues': [[24, 'P', 'H', [0, 4]], [25, 'E', 'H', [2, 4]]}
    :param ss_id: helix_1
    :param protein_name: AF-A0A009EPY2-F1-model_v4_TED01
    :param fig_dir: directory to save figure
    :return: HTML file
    """

    res_ids = [res[0] for res in group]
    residues = [res[1] for res in group]

    # move points to zero
    mean_coords = np.mean(points, axis=0)
    centered_coords = points - mean_coords

    # using PCA to calculate the vector
    pca = PCA(n_components=1)  # n_components=1 means line fitting
    pca.fit(centered_coords)

    # regression line direction
    line_direction = pca.components_[0]

    line_start = np.mean(centered_coords, axis=0)

    # calculate the distance between the points and the start point of the line
    distances = np.linalg.norm(centered_coords - line_start, axis=1)
    max_distance_idx = np.argmax(distances)

    # generate the two endpoints of the regression line
    t_min = -distances[max_distance_idx]
    t_max = distances[max_distance_idx]

    # calculate the two endpoints of the regression line
    line_end_1 = line_start + t_min * line_direction
    line_end_2 = line_start + t_max * line_direction

    # calculate the length of the regression line
    line_length = np.linalg.norm(line_end_2 - line_end_1)

    # modify the direction vector to the length of the regression line
    normalized_direction = line_direction / np.linalg.norm(line_direction)  # 单位化方向向量
    direction_vector = normalized_direction * line_length / 2  # 调整方向向量的长度为回归直线的长度

    # create a plotly figure
    fig = go.Figure()

    # draw
    fig.add_trace(go.Scatter3d(
        x=centered_coords[:, 0], y=centered_coords[:, 1], z=centered_coords[:, 2],
        mode='markers+text', marker=dict(size=5, color='blue'),
        text=[f'{res_ids[i]}{residues[i]}' for i in range(len(res_ids))],
        textposition='top center', name='Points'
    ))

    # draw line
    fig.add_trace(go.Scatter3d(
        x=[line_end_1[0], line_end_2[0]], y=[line_end_1[1], line_end_2[1]], z=[line_end_1[2], line_end_2[2]],
        mode='lines', line=dict(color='red', width=2),
        name='Regression Line'
    ))

    # draw direction
    fig.add_trace(go.Cone(
        x=[line_start[0]], y=[line_start[1]], z=[line_start[2]],
        u=[direction_vector[0]], v=[direction_vector[1]], w=[direction_vector[2]],
        colorscale='Viridis', showscale=False, sizemode="scaled", sizeref=0.2
    ))

    # calc range of data
    axis_min = np.array([centered_coords[:, 0].min(), centered_coords[:, 1].min(), centered_coords[:, 2].min()]).min()
    axis_max = np.array([centered_coords[:, 0].max(), centered_coords[:, 1].max(), centered_coords[:, 2].max()]).max()
    x_min, y_min, z_min = axis_min, axis_min, axis_min
    x_max, y_max, z_max = axis_max, axis_max, axis_max

    # layout
    margin = 1.1

    # setup labels
    fig.update_layout(
        title=f'{protein_name}: {ss_id}',
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            xaxis=dict(
                range=[x_min - margin, x_max + margin],
                showgrid=True,
                zeroline=False,
                showline=True,
                linewidth=2,
                linecolor='black',
                showticklabels=True,
                tickmode='auto',
                ticks='outside',
                ticklen=5,
            ),
            yaxis=dict(
                range=[y_min - margin, y_max + margin],
                showgrid=True,
                zeroline=False,
                showline=True,
                linewidth=2,
                linecolor='black',
                showticklabels=True,
                tickmode='auto',
                ticks='outside',
                ticklen=5,
            ),
            zaxis=dict(
                range=[z_min - margin, z_max + margin],
                showgrid=True,
                zeroline=False,
                showline=True,
                linewidth=2,
                linecolor='black',
                showticklabels=True,
                tickmode='auto',
                ticks='outside',
                ticklen=5,
            ),
            aspectmode="cube"
        ),
        margin=dict(l=10, r=10, b=10, t=40),  # 设置边距
        title_x=0.5,
    )

    # fig.show()
    fig_path = os.path.join(fig_dir, f'{protein_name}_{ss_id}.html')
    fig.write_html(fig_path)
    logger.info(f'Saved figure for {protein_name}_{ss_id} to dir: {fig_dir}')


def plot_vectors(vector1, vector2):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(0, 0, 0, color="black", s=50)  # 原点

    ax.quiver(0, 0, 0, *vector1, color='r', label="Vector 1", linewidth=2)
    ax.quiver(0, 0, 0, *vector2, color='b', label="Vector 2", linewidth=2)

    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    ax.legend()
    plt.show()
