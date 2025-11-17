# """
# # ------------------------------------------------------------------------------
# # Author:    Dou Zhixin
# # Email:     bj600800@gmail.com
# # DATE:      2023/07/19
#
# # Description: Service for hbond analysis.
# Pdb2Pqr for add hydrogen atom, biotite_hbond for calculate h-bonds
# # ------------------------------------------------------------------------------
# """
# import os
# import subprocess
# import tempfile
# import numpy as np
#
# import biotite.structure as struc
# import biotite.structure.io as strucio
# from biotite.structure import hbond, angle, distance
# import warnings
# warnings.filterwarnings("ignore")
#
#
# def add_h(input_pdb):
#     pqr_temp_file = tempfile.NamedTemporaryFile(delete=False)
#     pqr_temp_file_path = pqr_temp_file.name
#
#     pdb_addh_temp_file = tempfile.NamedTemporaryFile(delete=False)
#     pdb_addh_temp_file_path = pdb_addh_temp_file.name + ".pdb"
#
#     pdb2pqr_cmd = ["pdb2pqr", "--ff=AMBER", "--log-level", "CRITICAL", "--titration-state-method", "propka",
#                    "--pdb-output", pdb_addh_temp_file_path, input_pdb, pqr_temp_file_path]
#     pqr_temp_file.close()
#     os.remove(pqr_temp_file_path)
#     subprocess.run(pdb2pqr_cmd, check=True)
#     return pdb_addh_temp_file_path
#
#
# def calc_hbond(pdb, pdb_h):
#     stack = strucio.load_structure(pdb)
#     res_id, res_name = struc.get_residues(stack)
#     length = len(res_id)
#     stack_h = strucio.load_structure(pdb_h)
#     triplets_h = hbond(stack_h)
#     donor_acceptor = [(list(stack_h[bond].res_id[::2]), list(stack_h[bond].atom_name[::2])) for bond in triplets_h]
#     triplets = []
#     for resID_atmN in donor_acceptor:
#         res_id = resID_atmN[0]
#         atom_name = resID_atmN[1]
#         pair_idx = []
#         for i in range(0, 2):
#             hbond_mask = (stack.res_id == int(res_id[i])) & (stack.atom_name == atom_name[i])
#             pair_idx.append(np.where(hbond_mask)[0])
#         pair_idx = np.concatenate(pair_idx)
#         triplets.append(pair_idx)
#     triplets = np.array(triplets)
#     donor = stack_h[triplets_h[:, 0]]
#     donor_h = stack_h[triplets_h[:, 1]]
#     acceptor = stack_h[triplets_h[:, 2]]
#     theta = angle(donor, donor_h, acceptor)
#     degree = np.rad2deg(theta)
#     dist = distance(donor_h, acceptor)
#     zip_results = [*zip(degree, dist, triplets)]
#     hbond_results = len([tuple(i) for i in zip_results])/length
#     # for i in hbond_results:
#     #     print('+'.join([str(j) for j in list(stack[i[2]].res_id)]))
#     return hbond_results
#
#
# def analyze(structure_path):
#     pdb_addh_temp_file_path = add_h(structure_path)
#     hbond_results = calc_hbond(structure_path, pdb_addh_temp_file_path)
#     os.remove(pdb_addh_temp_file_path)
#     return hbond_results

import os
import subprocess
import tempfile
import warnings

import numpy as np
import biotite.structure as struc
import biotite.structure.io as strucio
from biotite.structure import hbond, angle, distance

warnings.filterwarnings("ignore")


# ------------------------------------------------------------------------------
# Step 1 — Add hydrogens via pdb2pqr
# ------------------------------------------------------------------------------
def add_hydrogens(input_pdb):
    """
    Use pdb2pqr to protonate the structure and output a temporary pdb file.
    """

    # Temporary files
    pqr_temp = tempfile.NamedTemporaryFile(delete=False)
    pqr_temp_path = pqr_temp.name
    pqr_temp.close()

    pdb_out_temp = tempfile.NamedTemporaryFile(delete=False)
    pdb_out_path = pdb_out_temp.name + ".pdb"
    pdb_out_temp.close()

    # Using updated options:
    # --nodebump: avoid error when structure has gaps or incomplete residues
    # --protonate-all: ensure all H atoms are added
    pdb2pqr_cmd = [
        "pdb2pqr",
        "--nodebump",
        "--protonate-all",
        "--ff=AMBER",
        "--log-level", "CRITICAL",
        "--titration-state-method", "propka",
        "--pdb-output", pdb_out_path,
        input_pdb,
        pqr_temp_path
    ]

    # Execute
    subprocess.run(pdb2pqr_cmd, check=True)

    # Remove the temporary pqr file
    os.remove(pqr_temp_path)

    return pdb_out_path


# ------------------------------------------------------------------------------
# Step 2 — Calculate H-bonds
# Unified logic combining both old and new behavior
# ------------------------------------------------------------------------------
def calculate_hbonds(pdb_original, pdb_h, return_mode="frequency"):
    """
    Unified hydrogen bond calculation:

    return_mode:
        "frequency" → return (#hbonds / length)
        "pairs"     → return list of donor/acceptor atom/residue pairs
    """

    # Original structure (without hydrogens)
    stack = strucio.load_structure(pdb_original)
    residue_ids, _ = struc.get_residues(stack)
    length = len(residue_ids)

    # Hydrogen-added structure
    stack_h = strucio.load_structure(pdb_h)

    # Biotite hydrogen bond detection
    triplets_h = hbond(stack_h)  # triplet: (donor_idx, H_idx, acceptor_idx)

    if return_mode == "pairs":

        hbonds = []
        # Format: {(atom_id, res_id), (atom_id, res_id)}

        for bond in triplets_h:
            # donor / acceptor atoms in the H-added structure
            donor = bond[0]
            acceptor = bond[2]

            donor_res = int(stack_h[donor].res_id)
            acceptor_res = int(stack_h[acceptor].res_id)

            donor_atom = donor + 1         # +1 to match user version
            acceptor_atom = acceptor + 1

            pair = {(donor_atom, donor_res), (acceptor_atom, acceptor_res)}

            if pair not in hbonds:
                hbonds.append(list(pair))

        return hbonds

    # --------------------------------------------------------------------------
    # Old behavior: compute frequency using original atom indices
    # --------------------------------------------------------------------------
    elif return_mode == "frequency":
        donor_acceptor = [(list(stack_h[bond].res_id[::2]),
                           list(stack_h[bond].atom_name[::2])) for bond in triplets_h]

        triplets = []
        for res_id_list, atom_name_list in donor_acceptor:
            idx_pair = []
            for i in range(2):
                mask = (stack.res_id == int(res_id_list[i])) & \
                       (stack.atom_name == atom_name_list[i])
                idx_pair.append(np.where(mask)[0])
            idx_pair = np.concatenate(idx_pair)
            triplets.append(idx_pair)

        triplets = np.array(triplets)
        donor = stack_h[triplets_h[:, 0]]
        donor_h = stack_h[triplets_h[:, 1]]
        acceptor = stack_h[triplets_h[:, 2]]

        theta = angle(donor, donor_h, acceptor)
        degree = np.rad2deg(theta)
        dist = distance(donor_h, acceptor)

        zip_results = list(zip(degree, dist, triplets))

        # frequency = #hbonds / protein_length
        return len(zip_results) / length

    else:
        raise ValueError("return_mode must be 'frequency' or 'pairs'")


# ------------------------------------------------------------------------------
# Step 3 — Unified high-level interface
# ------------------------------------------------------------------------------
def run(structure_path, return_mode="frequency"):
    """
    High-level unified H-bond analysis function.
    Adds hydrogens, computes H-bonds, then removes temporary files.
    """

    pdb_h = add_hydrogens(structure_path)
    try:
        result = calculate_hbonds(structure_path, pdb_h, return_mode=return_mode)
    finally:
        # Always remove hydrogen-added temporary pdb
        if os.path.exists(pdb_h):
            os.remove(pdb_h)

    return result

