"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/09/19

# Description: 
# ------------------------------------------------------------------------------
"""
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances
from scipy.stats import ttest_ind
import mplcursors
from sklearn.manifold import TSNE
from biotite.structure import gyration_radius, distance
from qprotein.utilities import align_map


def plot_heatmap(distance_matrix):
    plt.imshow(distance_matrix, cmap='RdYlBu', interpolation='nearest')
    plt.colorbar(label='Distance')
    plt.xlabel('Atom Index')
    plt.ylabel('Atom Index')
    plt.title('Distance Matrix Heatmap')
    plt.show()


def get_distance_matrix(atom_array, struct_id_map, max_id=363):
    ca_array = atom_array[(atom_array.atom_name == "CA")]
    distance_matrix = np.zeros((max_id, max_id))
    for i in range(len(ca_array)):
        for j in range(len(ca_array)):
            if i != j:
                distance_matrix[struct_id_map[i+1]-1][struct_id_map[j+1]-1] = distance(ca_array[i], ca_array[j])
    normalized_matrix = np.linalg.norm(distance_matrix)
    distance_matrix_normalized = distance_matrix / normalized_matrix
    return distance_matrix_normalized


def get_one_hot_cluster(hydrophobic_clusters, struct_id_map, max_id=363):
    one_hot_cluster = []
    for ori_cluster in hydrophobic_clusters:
        aligned_cluster = np.zeros((max_id,), dtype=int)
        map_cluster = list(map(lambda x: struct_id_map.get(x), ori_cluster))
        aligned_cluster[np.array(map_cluster) - 1] = 1
        one_hot_cluster.append([aligned_cluster, ori_cluster])
    return one_hot_cluster


def feature_cluster(atom_array, struct_name, hydrophobic_clusters, struct_id_map, max_id):
    one_hot_cluster = get_one_hot_cluster(hydrophobic_clusters, struct_id_map, max_id)
    distance_matrix = get_distance_matrix(atom_array, struct_id_map, max_id)
    feature = []
    for cluster in one_hot_cluster:
        vector = np.dot(cluster[0], distance_matrix)
        feature.append([[struct_name, cluster[1]], vector])
    return feature
        

def kmeans_cluster(all_feature, k):
    data = [one_cluster for example_dir in all_feature for one_structure in example_dir for one_cluster in one_structure]
    names = []
    vectors = []
    for cluster in data:
        names.append(cluster[0])
        vectors.append(cluster[1])
    vectors = np.array(vectors)
    # vectors = vectors / np.linalg.norm(vectors)
    kmeans = KMeans(n_clusters=k, init='k-means++')
    kmeans.fit(vectors)
    labels = kmeans.labels_
    # 输出聚类结果
    clusters = [[] for _ in range(k)]
    for i, label in enumerate(labels):
        clusters[label].append(names[i])
    tsne = TSNE(n_components=2, random_state=42, learning_rate=0.0001)
    X_tsne = tsne.fit_transform(np.array(vectors))

    check = []
    for i, x in enumerate(X_tsne):
        check.append([int(labels[i]), names[i][0], [str(j) for j in names[i][1]]])
    sorted_list = sorted(check, key=lambda x: x[0])
    for i in sorted_list:
        print(i[:2])
        print('select res ' + '+'.join(i[2]))

    return X_tsne, labels, names, clusters

def draw_cluster_plot(X_tsne, labels, names, clusters, k):
    fig, ax = plt.subplots()

    colors = ['C{}'.format(i) for i in range(k)]  # 颜色列表
    # 绘制散点图和标签
    for i in range(len(X_tsne)):
        x = X_tsne[i][0]
        y = X_tsne[i][1]
        label = labels[i]
        name = names[i]
        color = colors[label]
        ax.scatter(x, y, color=color, label=name, s=30)

    # 创建mplcursors对象
    mpl_cursors = mplcursors.cursor(hover=True)

    # 配置鼠标悬停时的标签显示
    @mpl_cursors.connect("add")
    def on_add(sel):
        label = sel.artist.get_label()
        sel.annotation.set_text(label)
        sel.annotation.get_bbox_patch().set(fc="white")


    # 创建图例
    color_patches = [mpatches.Patch(color=colors[i], label='Cluster '+str(i+1)) for i, cluster in enumerate(clusters)]
    fig.set_size_inches(8, 6)
    ax.legend(handles=color_patches, loc='upper right')
    # plt.savefig(r"D:\subject\active\1-qProtein\data\figure\result\gh11_hydrophobic.eps", format='eps', bbox_inches='tight')
    plt.show()

import os
from tqdm import tqdm
import biotite.structure.io as strucio
from qprotein.characterization.hydrophobic import analyze
structure_dir = [r"D:\subject\active\1-qProtein\data\enzymes\GH12\2_StrucMapping\positive",
                 r"D:\subject\active\1-qProtein\data\enzymes\GH12\2_StrucMapping\negative"]

all_ret = []
max_hydrophobic_cluster_lst = []
struct_id_map = pairwise.run()
max_id = max([v for key, value in struct_id_map.items() for k,v in value.items()])
for sdir in structure_dir:
    one_ret = []
    for struct in tqdm(os.listdir(sdir)):
        structure_path = os.path.join(sdir, struct)
        structure_name = os.path.splitext(os.path.basename(structure_path))[0]
        id_map = struct_id_map[structure_name]
        structure = strucio.load_structure(structure_path)
        out = analyze(atom_array=structure)
        hydrophobic_clusters = out['cluster_res']
        max_cluster = max(hydrophobic_clusters, key=len)
        max_hydrophobic_cluster_lst.append([structure_name, max_cluster])
        feature = feature_cluster(atom_array=structure, struct_name=structure_name, hydrophobic_clusters=hydrophobic_clusters, struct_id_map=id_map, max_id=max_id)
        one_ret.append(feature)
    all_ret.append(one_ret)


X_tsne, labels, names, clusters = kmeans_cluster(all_feature=all_ret, k=3)
draw_cluster_plot(X_tsne=X_tsne, labels=labels, names=names, clusters=clusters, k=3)
