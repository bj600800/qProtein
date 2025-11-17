# """
# # ------------------------------------------------------------------------------
# # Author:    Dou Zhixin
# # Email:     bj600800@gmail.com
# # DATE:      2024/11/27
#
# # Description: disulfide bond calculation and return the frequency of disulfide bonds.
# Original code source by Sphinx-Gallery: https://www.biotite-python.org/examples/gallery/structure/disulfide_bonds.html
# # ------------------------------------------------------------------------------
# """
# import biotite.structure as struc
# import numpy as np
#
#
# def run(structure, distance=2.05, distance_tol=3, dihedral=90, dihedral_tol=15):
#     """
#     largest distance for the abnormal ones
#     cite: https://doi.org/10.1093/protein/13.10.679
#     """
#     disulfide_bonds = []
#     sulfide_mask = (structure.res_name == "CYS") & \
#                    (structure.atom_name == "SG")
#     cell_list = struc.CellList(
#         structure,
#         cell_size=distance + distance_tol,
#         selection=sulfide_mask
#     )
#     for sulfide_i in np.where(sulfide_mask)[0]:
#         potential_bond_partner_indices = cell_list.get_atoms(
#             coord=structure.coord[sulfide_i],
#             radius=distance + distance_tol
#         )
#         for sulfide_j in potential_bond_partner_indices:
#             if sulfide_i == sulfide_j:
#                 continue
#             sg1 = structure[sulfide_i]
#             sg2 = structure[sulfide_j]
#             cb1 = structure[
#                 (structure.chain_id == sg1.chain_id) &
#                 (structure.res_id == sg1.res_id) &
#                 (structure.atom_name == "CB")
#                 ]
#             cb2 = structure[
#                 (structure.chain_id == sg2.chain_id) &
#                 (structure.res_id == sg2.res_id) &
#                 (structure.atom_name == "CB")
#                 ]
#             bond_dist = struc.distance(sg1, sg2)
#             bond_dihed = np.abs(np.rad2deg(struc.dihedral(cb1, sg1, sg2, cb2)))
#             if (distance - distance_tol < bond_dist < distance + distance_tol and
#                     dihedral - dihedral_tol < bond_dihed < dihedral + dihedral_tol):
#                 bond_tuple = {(int(sulfide_i) + 1, int(structure[sulfide_i].res_id)),
#                               (int(sulfide_j) + 1, int(structure[sulfide_j].res_id))}
#                 if bond_tuple not in disulfide_bonds:
#
#                     disulfide_bonds.append(bond_tuple)
#     return disulfide_bonds
#
#
# if __name__ == '__main__':
#     import os
#     import biotite.structure.io as strucio
#     structure_path = r"/Users/douzhixin/Developer/qProtein2/Data/pdb/1yna.pdb"
#     stack_h = strucio.load_structure(structure_path)
#     ret = run(stack_h)
#     print(ret)

# ------------------------------------------------------------------------------
# Author:      Dou Zhixin
# Email:       bj600800@gmail.com
# DATE:        2025/01/01
#
# Description:
#   Unified disulfide bond detection module.
#   Supports:
#       (1) returning disulfide bond frequency (old version)
#       (2) returning all disulfide bond pairs (new version)
#
# Original geometric criteria reference:
#   Sphinx-Gallery (Biotite)
#   https://www.biotite-python.org/examples/gallery/structure/disulfide_bonds.html
#
# Extended distance tolerance reference:
#   https://doi.org/10.1093/protein/13.10.679
# ------------------------------------------------------------------------------

import numpy as np
import biotite.structure as struc


def find_disulfide_bonds(
        structure,
        distance=2.05,
        distance_tol=0.05,
        dihedral=90,
        dihedral_tol=15,
        return_mode="frequency",
        extended_tol=False
):
    """
    Detect disulfide bonds in a protein structure.

    Parameters
    ----------
    structure : AtomArray or AtomArrayStack
        Biotite structure object.
    distance : float
        Ideal SG–SG bond distance.
    distance_tol : float
        Allowed deviation from ideal bond distance.
    dihedral : float
        Expected dihedral angle CB–SG–SG–CB.
    dihedral_tol : float
        Allowed deviation from ideal dihedral angle.
    return_mode : {"frequency", "bonds"}
        - "frequency": return number_of_bonds / number_of_residues
        - "bonds":     return list of disulfide bond tuples
    extended_tol : bool
        Whether to use extended tolerance (±3 Å) for abnormal bonds
        based on https://doi.org/10.1093/protein/13.10.679.

    Returns
    -------
    float or list
        Frequency or list of disulfide bond index/residue pairs.
    """

    # If user chooses extended tolerance
    if extended_tol:
        distance_tol = 3.0

    disulfide_bonds = []

    # Mask SG atoms in cysteine
    sulfide_mask = (structure.res_name == "CYS") & (structure.atom_name == "SG")

    # Speed up neighbor search
    cell_list = struc.CellList(
        structure,
        cell_size=distance + distance_tol,
        selection=sulfide_mask
    )

    sulfide_indices = np.where(sulfide_mask)[0]

    for i in sulfide_indices:
        partners = cell_list.get_atoms(
            coord=structure.coord[i],
            radius=distance + distance_tol
        )

        for j in partners:
            if i == j:
                continue  # skip itself

            sg1 = structure[i]
            sg2 = structure[j]

            # Retrieve CB atoms
            cb1 = structure[
                (structure.chain_id == sg1.chain_id) &
                (structure.res_id  == sg1.res_id) &
                (structure.atom_name == "CB")
            ]
            cb2 = structure[
                (structure.chain_id == sg2.chain_id) &
                (structure.res_id  == sg2.res_id) &
                (structure.atom_name == "CB")
            ]

            if len(cb1) == 0 or len(cb2) == 0:
                continue  # skip incomplete residues

            bond_dist = struc.distance(sg1, sg2)
            bond_dihed = abs(np.rad2deg(struc.dihedral(cb1, sg1, sg2, cb2)))

            # geometric criteria
            cond_dist = (distance - distance_tol) < bond_dist < (distance + distance_tol)
            cond_dihed = (dihedral - dihedral_tol) < bond_dihed < (dihedral + dihedral_tol)

            if cond_dist and cond_dihed:
                # Always store sorted tuple to avoid duplicates
                bond_tuple = tuple(sorted([
                    (int(i) + 1, int(sg1.res_id)),
                    (int(j) + 1, int(sg2.res_id))
                ]))

                if bond_tuple not in disulfide_bonds:
                    disulfide_bonds.append(bond_tuple)

    # Return mode selector
    if return_mode == "frequency":
        unique_res_ids = len(set(structure.res_id))
        return len(disulfide_bonds) / unique_res_ids if unique_res_ids > 0 else 0.0

    elif return_mode == "bonds":
        return disulfide_bonds

    else:
        raise ValueError("return_mode should be 'frequency' or 'bonds'")


# ------------------------------------------------------------------------------
# Main test
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    import biotite.structure.io as strucio

    structure_path = r"/Users/douzhixin/Developer/qProtein2/Data/pdb/1yna.pdb"
    structure = strucio.load_structure(structure_path)

    print("=== Disulfide bonds ===")
    ret_bonds = find_disulfide_bonds(structure, return_mode="bonds", extended_tol=True)
    print(ret_bonds)

    print("\n=== Frequency ===")
    ret_freq = find_disulfide_bonds(structure, return_mode="frequency")
    print(ret_freq)
