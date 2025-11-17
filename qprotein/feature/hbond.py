import os
import subprocess
import tempfile
import warnings

import numpy as np
import biotite.structure as struc
import biotite.structure.io as strucio
from biotite.structure import hbond, angle, distance

warnings.filterwarnings("ignore")

def add_hydrogens(input_pdb):
    """
    Use pdb2pqr to protonate the structure and output a temporary pdb file.
    """
    pqr_temp = tempfile.NamedTemporaryFile(delete=False)
    pqr_temp_path = pqr_temp.name
    pqr_temp.close()

    pdb_out_temp = tempfile.NamedTemporaryFile(delete=False)
    pdb_out_path = pdb_out_temp.name + ".pdb"
    pdb_out_temp.close()

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

    subprocess.run(pdb2pqr_cmd, check=True)

    os.remove(pqr_temp_path)

    return pdb_out_path

def calculate_hbonds(pdb_original, pdb_h, return_mode):
    """
    Unified hydrogen bond calculation:

    return_mode:
        "frequency" → return (#hbonds / length)
        "pairs"     → return list of donor/acceptor atom/residue pairs
    """
    stack = strucio.load_structure(pdb_original)
    residue_ids, _ = struc.get_residues(stack)
    length = len(residue_ids)

    stack_h = strucio.load_structure(pdb_h)

    triplets_h = hbond(stack_h)  # triplet: (donor_idx, H_idx, acceptor_idx)

    if return_mode == "pairs":

        hbonds = []

        for bond in triplets_h:
            donor = bond[0]
            acceptor = bond[2]

            donor_res = int(stack_h[donor].res_id)
            acceptor_res = int(stack_h[acceptor].res_id)

            donor_atom = donor + 1
            acceptor_atom = acceptor + 1

            pair = {(donor_atom, donor_res), (acceptor_atom, acceptor_res)}

            if pair not in hbonds:
                hbonds.append(list(pair))

        return hbonds

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

        return len(zip_results) / length, length

    else:
        raise ValueError("return_mode must be 'frequency' or 'pairs'")

def run(structure_path, return_mode):
    """
    High-level unified H-bond analysis function.
    Adds hydrogens, computes H-bonds, then removes temporary files.
    """

    pdb_h = add_hydrogens(structure_path)
    try:
        result = calculate_hbonds(structure_path, pdb_h, return_mode=return_mode)
    finally:
        if os.path.exists(pdb_h):
            os.remove(pdb_h)

    return result

