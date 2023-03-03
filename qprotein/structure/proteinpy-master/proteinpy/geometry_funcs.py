## !/usr/bin/python3

from __future__ import print_function
# import pandas as pd
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R


def calc_angle(vec1, vec2):
    """Calculate the angle between two vectors."""
    if not (int(np.linalg.norm(vec1)*np.linalg.norm(vec2))>0):
        return -1
    cos_val = np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
    if not (cos_val>=-1.0 and cos_val<=1.0):
        return -1
    angle = np.arccos(np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))
    return np.degrees(angle)

def calc_dihedral_angle(vec1, vec2, vec3):
    """Calculate the dihedral angle between two planes formed by three vectors plane(vec1-vec2) ~ plane(vec1-vec3)."""
    pln_vec1 = np.cross(vec1, vec2)
    pln_vec2 = np.cross(vec1, vec3)
    dihedral_angle = calc_angle(pln_vec1, pln_vec2)
    return dihedral_angle


def calc_tetrahedron_angles(vecs):
    """
    Returns a list of angles between any two vectors it takes in.
    """
    angles = []
    for i in range(len(vecs)):
        for j in range(i+1, len(vecs)):
            angles.append(calc_angle(vecs[i], vecs[j]))
    
    return angles

def calc_tetrahedron_dihedral_angles(vecs):
    """
    Returns a list of dihedral angles between any three vectors it takes in.
    """
    dihedral_angles = []
    for i in range(len(vecs)):
        for j in range(len(vecs)):
            for k in range(len(vecs)):
                if i!=j and i!=k and j<k:
                    dihedral_angles.append(calc_dihedral_angle(vecs[i], vecs[j], vecs[k]))
    
    return dihedral_angles


def calc_cycle_angles(vecs):
    """
    Returns a list of obtuse angles between any two neighbor vectors from a list of cycled vectors it takes in
    """
    angles = []
    for i in range(len(vecs)-1):
        angles.append(calc_angle(vecs[i], -vecs[i+1]))
    angles.append(calc_angle(vecs[0], -vecs[-1]))

    return angles

def calc_cycle_dihedral_angles(vecs):
    """
    Returns a list of dihedral angles between any two neighbor vectors from a list of cycled vectors it takes in
    """
    cycle_dihedral_angles = []
    for i in range(len(vecs)-1):
        cycle_dihedral_angles.append(calc_dihedral_angle(-vecs[i], vecs[i+1], vecs[i-1]))
    cycle_dihedral_angles.append(calc_dihedral_angle(-vecs[len(vecs)-1], vecs[0], vecs[len(vecs)-2]))

    return cycle_dihedral_angles

def draw_tetrahedron(fig, vecs):
    """plot a 3D tetrahedron using given 4 vectors."""

    ax = fig.gca(projection='3d')
    colors = ['red', 'green', 'blue', 'purple']
    for i, vec in enumerate(vecs):
        x_line = np.linspace(0, vec[0], 100)
        y_line = np.linspace(0, vec[1], 100)
        z_line = np.linspace(0, vec[2], 100)
        ax.plot3D(x_line, y_line, z_line, color=colors[0])
    # plt.show()

def match_tetrahedron(vecs1, vecs2):
    """
    Rotate a set of 4 vectors (vecs2) to match another set of 4 vectors (vecs1).
    Returns: rotated vecs2.
    """
    rotation, __ = R.match_vectors(vecs1, vecs2)
    return rotation.apply(vecs2)

def superimpose_tetrahedrons(listofvecs):
    """superimpose a list of different tetrahedrons (in the form of 4 vectors), and plot 3D figure of the superimposition."""
    
    fig = plt.figure()
    vecs0 = np.mean(np.array(listofvecs), axis=0)#np.mean(listofvecs)

    for vecs in listofvecs:
        rotated_vecs = match_tetrahedron(vecs0, vecs)
        draw_tetrahedron(fig, rotated_vecs)
    plt.show()

# v1 = np.array([1, 0, 0])
# v2 = np.array([0, 1, 0])
# v3 = np.array([0, 0, 1])
# v4 = np.array([-1, -1, -1])

# print(calc_dihedral_angles([v1, v2, v3, v4]))

# v1 = np.array([-0.5, 0.5, 0.5])
# v2 = np.array([0.5, 0.5, -0.5])
# v3 = np.array([1, 0, 0])
# v4 = np.array([0.5, -0.5, -0.5])
# v5 = np.array([-0.5, -0.5, 0.5])
# v6 = np.array([-1, 0, 0])
# print(calc_cycle_dihedral_angles([v1, v2, v3, v4, v5, v6]))

# p1 = [0, 0, 0]

# p2 = np.array([1, 0, 0])
# p3 = np.array([0, 1, 0])
# p4 = np.array([0, 0, 1])
# p5 = np.array([-0.577, -0.577, -0.577])
# vs1 = np.array([p2,p3,p4,p5])

# p2_ = np.array([0.707, -0.707, 0])
# p3_ = np.array([0.707, 0.707, 0])
# p4_ = np.array([0, 0, 1])
# p5_ = np.array([-0.817, 0, -0.577])
# vs2 = np.array([p2_,p3_,p4_,p5_])

# vs3 = np.array([[np.sqrt(8.0/9.0), 0, -1.0/3.0], \
#                 [-np.sqrt(2.0/9.0), np.sqrt(2.0/3.0), -1.0/3.0], \
#                 [-np.sqrt(2.0/9.0), -np.sqrt(2.0/3.0), -1.0/3.0], \
#                 [0, 0, 1]])

# # draw_tetrahedron(p2, p3, p4, p5)
# superimpose_tetrahedrons(vs1,)