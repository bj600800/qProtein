from tempfile import gettempdir
import biotite
import biotite.structure as struc
import biotite.structure.io as strucio
import biotite.database.rcsb as rcsb
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


file_name = rcsb.fetch("1aki", "mmtf", gettempdir())
array = strucio.load_structure(file_name)
# We only consider CA atoms
ca = array[array.atom_name == "CA"]
# 7 Angstrom adjacency threshold
threshold = 7
# Create cell list of the CA atom array
# for efficient measurement of adjacency
cell_list = struc.CellList(ca, cell_size=threshold)
adjacency_matrix = cell_list.create_adjacency_matrix(threshold)

figure = plt.figure()
ax = figure.add_subplot(111)
cmap = ListedColormap(["white", biotite.colors["dimgreen"]])
ax.matshow(adjacency_matrix, cmap=cmap, origin="lower")
ax.xaxis.tick_bottom()
ax.set_aspect("equal")
ax.set_xlabel("Residue number")
ax.set_ylabel("Residue number")
ax.set_title("Adjacency matrix of the lysozyme crystal structure")
figure.tight_layout()
plt.show()