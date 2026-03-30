"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2025/3/26

# Description: Patch program for get_pdb.py
Split protein domains from src and save target domain to dir
# ------------------------------------------------------------------------------
"""
import argparse
import os
import subprocess
from pathlib import Path

from tqdm import tqdm

from qpacking.utils import logger

logger = logger.setup_log(name=__name__)

#### ARGUMENTS PARSER ####
parser = argparse.ArgumentParser(description='domain segmentation')
parser.add_argument('--id', required=True, help='output id file')
parser.add_argument('--src', required=True, help='source structure dir')
parser.add_argument('--tar', required=True, help='target domain structure dir')
args = parser.parse_args()


#### END OF ARGUMENTS PARSER ####

def get_ted_ids(ted_domain_file):
    """
    Get protein ID and domain segmentation residue positions from ted_domain_file.

    :param ted_domain_file: Path to the file containing two columns: first protein, second domain segmentation.
    :return: Dictionary with protein_domain_id as keys and residue start-end positions as values.
    """
    if os.path.exists(ted_domain_file):
        with open(ted_domain_file, "r") as f:
            id_dict = {line.split('\t')[0].strip(): line.split('\t')[1].strip() for line in f}

        logger.info(f"Found {len(id_dict)} items in {ted_domain_file}")
    else:
        logger.error(f"No id files")
        exit(1)
    return id_dict


def run(id_dict, src_dir, domain_dir):
    """
    Split domains from source structure files and save to target directory.

    :param id_dict: Dictionary with protein_domain_id as keys and residue start-end positions as values.
    :param src_dir: Path to the source structure dir.
    :param domain_dir: Path to the domain dir.
    """
    if not os.path.exists(src_dir):
        logger.info(f"src structure dir doesn't exist")
        return
    error_log = os.path.join(os.path.dirname(src_dir), 'domain_' + os.path.basename(src_dir) + '.log')

    os.makedirs(domain_dir, exist_ok=True)
    pdb_files = list(Path(src_dir).glob("*.pdb"))
    for src_structure in tqdm(pdb_files, desc='Splitting domains'):
        protein_name = src_structure.stem
        domain_path = os.path.join(domain_dir, protein_name + ".pdb")

        if os.path.exists(domain_path):
            logger.info(f"Skipping {protein_name}, already processed.")
            continue

        res_pos = id_dict.get(protein_name)
        if not res_pos:
            logger.warning(f"No domain information for {protein_name}, skipping.")
            continue
        res_pos = res_pos.replace('-', ':')
        try:
            command = ['pdb_selres', f'-{res_pos}', str(src_structure)]
            result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            with open(domain_path, "w") as f:
                f.write(result.stdout.decode())
        except Exception as e:
            with open(error_log, "a") as f:
                f.write(protein_name + '\t' + res_pos + '\n')
            logger.error(f"Error processing {protein_name}: {e}")


if __name__ == '__main__':
    ted_domain_file = args.id
    src_structure_dir = args.src
    target_structure_dir = args.tar
    id_dict = get_ted_ids(ted_domain_file)
    run(id_dict, src_structure_dir, target_structure_dir)
