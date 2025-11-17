"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2025/1/7

# Description: multiple structure alignment followed with phylogenetic analysis
including site-specific conservation, evolution-function relationship.
Output results as landscape plot, pymol visualization and excel details.
figure name: conservative landscape and weblogo. row name: site, column name: aa
use grid, left landscape right weblogo，up label count for
# ------------------------------------------------------------------------------
"""
import os
import subprocess
import time

import biotite.structure as struc
import biotite.structure.io as strucio
import pandas as pd

from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


def read_labels(label_file):
    label_dict = {}
    with open(label_file, 'r') as f:
        for line in f:
            split_parts = line.rstrip().split()
            if len(split_parts) == 2:
                protein_name, label = split_parts
                label_dict[protein_name] = float(label)
    return label_dict


def align_pdbs(pdb_dir, template_name, label_dict, usalign_bin):
    """
    Align all structures

    Args:
        pdb_dir: pdb dir path
        usalign_bin: usalign binary path

    Returns: Path of aligned results as a fasta file.

    """
    time1 = time.time()
    logger.info("Aligning PDBs...")

    chain_list = [file for file in os.listdir(pdb_dir) if file.endswith(".pdb")]
    total_protein_num = len(chain_list)

    chain_list_path = os.path.join(pdb_dir, 'chain_list')
    with open(chain_list_path, 'w') as f:
        f.writelines(f"{line.rstrip('.pdb')}\n" for line in chain_list if line != template_name)

    template_pdb = os.path.join(pdb_dir, template_name)
    structure = strucio.load_structure(template_pdb)
    res_idx = list(set(structure[struc.filter_canonical_amino_acids(structure)].res_id.tolist()))

    align_seq_df = pd.DataFrame(columns=res_idx)

    usalign_cmd = [
        usalign_bin,
        "-mol",
        "prot",
        template_pdb,
        "-dir2",
        pdb_dir,
        chain_list_path,
        # "-do",  #show paired residues distances
        "-outfmt",
        "1"
    ]

    try:
        result = subprocess.run(usalign_cmd, check=True, capture_output=True, text=True)
        align_output = [align_pair.strip().split('\n') for align_pair in result.stdout.split("$$$$")[:-1]]
        align_res_list = [[out[1], out[3],
                           label_dict[out[2].split('\t')[0].split('/')[-1].split('.pdb')[0]]]
                          for out in align_output]

        time2 = time.time()
        logger.info(f"Aligning PDBs... DONE in {round(time2-time1, 2)} seconds.")

        return align_res_list, total_protein_num, res_idx, align_seq_df

    except subprocess.CalledProcessError as e:
        logger.error(e)


def init_aa():
    return {
        'A': [],
        'C': [],
        'D': [],
        'E': [],
        'F': [],
        'G': [],
        'H': [],
        'I': [],
        'K': [],
        'L': [],
        'M': [],
        'N': [],
        'P': [],
        'Q': [],
        'R': [],
        'S': [],
        'T': [],
        'V': [],
        'W': [],
        'Y': []
    }


def map_positions_with_gap(seq_with_gap, seq_without_gap):
    """
    map two sequences with and without gaps to each other
    1. seq_with_gap: A--BC-DEF
    2. seq_without_gap: ABCDEF
    3. mapping: {0:0, 3:1, 4:2, 6:3, 7:4, 8:5}
    4. return mapping

    Args:
        seq_with_gap: sequence with gap
        seq_without_gap: sequence without gap

    Returns:
        list:
    """

    mapping = {}
    index_seq_without_gap = 0

    for index_seq_with_gap in range(len(seq_with_gap)):
        if seq_with_gap[index_seq_with_gap] != '-':
            if index_seq_without_gap < len(seq_without_gap):
                mapping[index_seq_with_gap] = index_seq_without_gap
                index_seq_without_gap += 1

    return mapping


def create_landscape(template_label, align_res_list, align_seq_df):
    landscape = []

    # process template
    template_without_gaps = align_res_list[0][0].replace('-','')
    align_seq_df.loc[0] = list(template_without_gaps)

    for res in template_without_gaps:
        aa_dict = init_aa()
        aa_dict[res].append(template_label)
        landscape.append(aa_dict)

    # continue to process candidates with landscape list
    for align_pair in align_res_list:
        candi_res_str = ''
        template, candidate = align_pair[:2]
        mapping = map_positions_with_gap(template, template_without_gaps)

        for i in range(len(template)):
            if template[i] != '-':
                pos_without_gap = mapping[i]
                candi_res = candidate[i]
                candi_res_str += candi_res
                candi_label = align_pair[2]
                if candi_res != '-':
                    landscape[pos_without_gap][candi_res].append(candi_label)
        align_seq_df.loc[len(align_seq_df)] = list(candi_res_str)

    return landscape, align_seq_df

def calc_conservation(landscape, total_protein_num, res_idx, output_xlsx, label_pos_thrshold):
    """
    Calculate residue conservation for all positions
    Args:
        landscape: raw data of protein evolution with labels
        total_protein_num: sample size
        res_idx: [1,2,3...100]

    Returns: df_landscape: dataframe with residue conservation

    """

    label_pos_thrshold = float(label_pos_thrshold)
    amino_acids = ['S', 'T', 'N', 'Q', 'R', 'K', 'H', 'E', 'D', 'P',
                   'G', 'A', 'V', 'I', 'L', 'C', 'M', 'F', 'Y', 'W']

    df = pd.DataFrame(index=amino_acids, columns=res_idx)

    for pos_idx in range(len(df.columns)):

        for aa in df.index:
            labels = landscape[pos_idx][aa]
            current_aa_conservation = round(len(landscape[pos_idx][aa]) / total_protein_num, 4)
            pos_label_ratio = round((pd.Series(labels) > label_pos_thrshold).mean(), 4)

            if current_aa_conservation:
                df.loc[aa, df.columns[pos_idx]] = [pos_label_ratio, current_aa_conservation]
            else:
                df.loc[aa, df.columns[pos_idx]] = None

    df = pd.DataFrame(index=amino_acids, columns=res_idx)

    for pos_idx in range(len(df.columns)):

        for aa in df.index:
            labels = landscape[pos_idx][aa]
            current_aa_conservation = round(len(landscape[pos_idx][aa]) / total_protein_num, 4)
            pos_label_ratio = round((pd.Series(labels) > label_pos_thrshold).mean(), 4)

            if current_aa_conservation:
                df.loc[aa, df.columns[pos_idx]] = [pos_label_ratio, current_aa_conservation]
            else:
                df.loc[aa, df.columns[pos_idx]] = None

    df_write = df.map(lambda x: f"{x[0]}, {x[1]}" if isinstance(x, list) else x)
    excel_file = output_xlsx
    df_write.to_excel(excel_file)

    return df


def run(label_file, pdb_dir, template_name, output_xlsx, label_pos_thrshold, usalign_bin):
    label_dict = read_labels(label_file)
    template_label = label_dict[template_name.rstrip('.pdb')]
    align_res_list, total_protein_num, res_idx, align_seq_df = align_pdbs(pdb_dir, template_name, label_dict, usalign_bin)
    landscape, aligned_seq_df = create_landscape(template_label, align_res_list, align_seq_df)
    df_landscape = calc_conservation(landscape, total_protein_num, res_idx, output_xlsx, label_pos_thrshold)
    return df_landscape, aligned_seq_df
