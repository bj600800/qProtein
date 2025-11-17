"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2024/11/25

# Description: Compute surface charged residues.
# ------------------------------------------------------------------------------
"""

import os
import tempfile
import warnings
from pathlib import Path
from subprocess import Popen, PIPE

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from tqdm import tqdm

from qpacking.common import logger

logger = logger.setup_log(name=__name__)
warnings.simplefilter('ignore', PDBConstructionWarning)


def get_dssp_dat(struct_path, dssp_bin):
    struct_file_name = os.path.split(struct_path)[1].split('.')[0]
    parser = PDBParser()
    model = parser.get_structure(struct_file_name, struct_path)[0]
    chain = next(model.get_chains())
    start_residue_num = next(chain.get_residues()).id[1]
    try:
        dssp = DSSP(model, struct_path, dssp=dssp_bin)
        dssp_dat = list(map(lambda residue_dat: (residue_dat[0] + start_residue_num - 1, *residue_dat[1:]), dssp))
        return dssp_dat
    except FileNotFoundError:
        logger.error('dssp not found. Please designate to the directory where dssp binary is located.')
        return None

def get_aa_charge(dssp_dat):
    residue_physical = {
        'positive': ('K', 'R', 'H'),
        'negative': ('D', 'E')
    }
    residues = [i[:4] for i in dssp_dat]
    total = len(residues)

    positive = 0
    negative = 0
    # check residue charge
    charge_position = []
    for res in residues:

        if res[1] != 'X' and res[3] > 0.05:
            if res[1] in residue_physical['positive']:
                positive += 1
                charge_position.append(res[:2])  # add positive charge
            elif res[1] in residue_physical['negative']:
                negative += 1
                charge_position.append(res[:2])  # add negative charge
    prop_positive = '{:.3}'.format(positive / total * 100)
    prop_negative = '{:.3}'.format(negative / total * 100)

    percent = {'negative': prop_negative, 'positive': prop_positive}
    return percent, charge_position


def get_charges(pdb_file, pdb2pqr_bin, apbs_bin):
    """
        Calls APBS, pdb2pqr
    """
    base_pdb = os.path.basename(pdb_file.stem)
    tmp_dir = tempfile.TemporaryDirectory().name
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    apbs_input = os.path.join(tmp_dir, base_pdb + '.in')
    pqr_file = os.path.join(tmp_dir, base_pdb + '.pqr')
    pdb2pqr_cmd = [
        pdb2pqr_bin,
        "--nodebump",
        "--ff=PARSE",
        "--protonate-all",
        "--whitespace",
        "--apbs-input",
        apbs_input,
        pdb_file,
        pqr_file

    ]

    try:
        p1 = Popen(pdb2pqr_cmd, stdout=PIPE, stderr=PIPE, cwd=tmp_dir)
        p1.communicate()

        apbs_cmd = [apbs_bin, os.path.join(tmp_dir, base_pdb + ".in")]
        p2 = Popen(apbs_cmd, stdout=PIPE, stderr=PIPE, cwd=tmp_dir)
        stdout, stderr = p2.communicate()
    except:
        raise

    if stderr:
        logger.error(stderr)
        return

    apbs_out = stdout.decode("utf-8").split("\n")
    charge = ''
    for line in apbs_out:
        if line.startswith("  Net charge"):
            charge = str(float(line.split(" ")[4]))

    return int(float(charge))


def run(structure_dir, dssp_bin, pdb2pqr_bin, apbs_bin):
    logger.info('Analyzing surface charged residues...')
    ret = []
    pdb_files = list(Path(structure_dir).glob("*.pdb"))
    for struc_path in tqdm(pdb_files, desc='Computing surface charged residues'):
        file_name = struc_path.stem
        dssp_dat = get_dssp_dat(struc_path, dssp_bin=dssp_bin)
        charge = get_charges(struc_path, pdb2pqr_bin=pdb2pqr_bin, apbs_bin=apbs_bin)
        if dssp_dat:
            result = get_aa_charge(dssp_dat)
            ret.append([file_name, result[1], result[0], charge])
    return ret

if __name__ == '__main__':
    structure_dir = r"/Users/douzhixin/Developer/qProtein2/Data/pdb"
    dssp_bin = "mkdssp"
    pdb2pqr_bin = "pdb2pqr"
    apbs_bin = "apbs"
    ret = run(structure_dir, dssp_bin, pdb2pqr_bin, apbs_bin)
    print(ret)
