"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/07/11

# Description: disulfide bond calculation and return the frequency of disulfide bonds.
Original code source by Sphinx-Gallery: https://www.biotite-python.org/examples/gallery/structure/disulfide_bonds.html
# ------------------------------------------------------------------------------
"""
import numpy as np
import biotite.structure as struc


def analyze(structure, distance=2.05, distance_tol=0.05, dihedral=90, dihedral_tol=15):
    disulfide_bonds = []
    length = len(set(structure.res_id))
    sulfide_mask = (structure.res_name == "CYS") & \
                   (structure.atom_name == "SG")
    cell_list = struc.CellList(
        structure,
        cell_size=distance + distance_tol,
        selection=sulfide_mask
    )
    for sulfide_i in np.where(sulfide_mask)[0]:
        potential_bond_partner_indices = cell_list.get_atoms(
            coord=structure.coord[sulfide_i],
            radius=distance + distance_tol
        )
        for sulfide_j in potential_bond_partner_indices:
            if sulfide_i == sulfide_j:
                continue
            sg1 = structure[sulfide_i]
            sg2 = structure[sulfide_j]
            cb1 = structure[
                (structure.chain_id == sg1.chain_id) &
                (structure.res_id == sg1.res_id) &
                (structure.atom_name == "CB")
                ]
            cb2 = structure[
                (structure.chain_id == sg2.chain_id) &
                (structure.res_id == sg2.res_id) &
                (structure.atom_name == "CB")
                ]
            bond_dist = struc.distance(sg1, sg2)
            bond_dihed = np.abs(np.rad2deg(struc.dihedral(cb1, sg1, sg2, cb2)))
            if (distance - distance_tol < bond_dist < distance + distance_tol and
                    dihedral - dihedral_tol < bond_dihed < dihedral + dihedral_tol):
                bond_tuple = tuple(sorted((sulfide_i, sulfide_j)))
                if bond_tuple not in disulfide_bonds:
                    disulfide_bonds.append(bond_tuple)
    return len(disulfide_bonds)/length


    