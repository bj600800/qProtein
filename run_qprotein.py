"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2025/01/01

# Description: qProtein run file
# ------------------------------------------------------------------------------
"""
import os
import argparse
import subprocess
import time
from tqdm import tqdm
import biotite.structure.io as strucio

from qprotein.seq2struct.sequence import get_id
from qprotein.seq2struct.afdb import get_struct
from qprotein.seq2struct import esmfold
from qprotein.feature import hydrophobic, hbond, salt_bridge, disulfide_bond
from qprotein.analysis import overall, local
from qprotein.utilities import logger

from qprotein.analysis import internal, surface, landscape
logger = logger.setup_log(name=__name__)

def crawl_struct(id_list, structure_folder):
    """Download AFDB structures unless already existing."""
    def search_exist_struct(structure_folder):
        exists = [file for file in os.listdir(structure_folder)
                  if os.path.isfile(os.path.join(structure_folder, file))
                  and os.path.getsize(os.path.join(structure_folder, file)) > 0]
        return [os.path.splitext(item)[0] for item in exists]

    id_list = list(set(id_list) - set(search_exist_struct(structure_folder)))
    if id_list:
        for uniprot_id in tqdm(id_list):
            pdb_string = get_struct(uniprot_id)
            if pdb_string:
                save_path = os.path.join(structure_folder, uniprot_id + ".pdb")
                with open(save_path, "w") as f:
                    f.write(pdb_string)

    logger.info(f"Crawled structures: {len(os.listdir(structure_folder))}")


def calc_feature(pdb_list, return_mode):
    logger.info(f"Calculating features, please wait...")
    feature_dict = {}
    for pdb_file in tqdm(pdb_list):
        name = os.path.basename(pdb_file).split(".")[0]
        structure = strucio.load_structure(pdb_file)
        hbond_freq, length = hbond.run(pdb_file, return_mode)
        salt_freq, length = salt_bridge.run(structure, return_mode)
        disul_freq, length = disulfide_bond.run(structure, return_mode)
        hydrophobic_ret = hydrophobic.run(structure)
        feature_dict[name] = {
            "hydrophobic": hydrophobic_ret,
            "hbond": hbond_freq,
            "saltbridge": salt_freq,
            "disulfide": disul_freq,
            "length": length
        }
    return feature_dict

def calc_local_hydrophobic(pdb_list):
    hydro_dict = {}
    for pdb_file in tqdm(pdb_list):
        name = os.path.basename(pdb_file).split(".")[0]
        structure = strucio.load_structure(pdb_file)
        hydrophobic_ret = hydrophobic.run(structure)
        hydro_dict[name] = hydrophobic_ret
    return hydro_dict

def align_structure(structure_folder, usalign_binary="usalign"):
    files = [file for file in os.listdir(structure_folder) if file.endswith(".pdb")]
    name_txt = os.path.join(structure_folder, "name.txt")

    with open(name_txt, "w") as f:
        for fn in files:
            f.write(fn + "\n")

    out_file = os.path.join(os.path.dirname(structure_folder), "usalign_out.fasta")

    cmd = [usalign_binary, "-dir", structure_folder, name_txt, "-suffix", ".pdb", "-mm", "4"]
    result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE)

    with open(out_file, "w") as f:
        f.write(result.stdout.decode())

    os.remove(name_txt)
    return out_file

def main():
    modes = ["fetch", "overall", "local", "visual", "surface", "landscape"]

    parser = argparse.ArgumentParser(description="qProtein analysis kit")
    parser.add_argument("--mode", type=str, choices=modes, required=True)
    parser.add_argument("--fasta")
    parser.add_argument("--id")
    parser.add_argument("--pre_pdb")
    parser.add_argument("--work_dir", required=True)
    parser.add_argument("--template_name")
    parser.add_argument("--dist1")
    parser.add_argument("--dist2")
    parser.add_argument("--config_file")
    parser.add_argument("--pml", action="store_true")
    parser.add_argument("--positions")
    parser.add_argument("--label_file")
    parser.add_argument("--label_pos_threshold")

    args = parser.parse_args()

    # ------------------------------------------------------------------------------
    # MODE — FETCH: Fetch pdbs from AFDB or ESMFold
    # ------------------------------------------------------------------------------
    if args.mode == "fetch":
        structure_dir = os.path.join(args.work_dir, "structure")
        os.makedirs(structure_dir, exist_ok=True)

        # 1. get structure
        if args.fasta:
            logger.info("Predicting structures by ESMFold")
            esmfold.run(
                fasta_file=args.fasta,
                structure_dir=structure_dir,
                esm_script="/opt/app/esm-main/scripts/fold.py",
                esm_dir="/opt/app/esm-main/",
            )

        elif args.id:
            logger.info("Crawling AlphaFold DB")
            crawl_struct(get_id(args.id), structure_dir)

        else:
            logger.error("fetch mode requires either --fasta or --id")
            logger.error("Usage example: --mode fetch --fasta xxx.fasta  OR  --mode fetch --id id_list.txt")
            return

        logger.info(f"Structure fetch completed. Saved to: {structure_dir}")
        return

    # ------------------------------------------------------------------------------
    # MODE — OVERALL: Calculate overall protein structure features
    # ------------------------------------------------------------------------------
    elif args.mode == "overall":
        start = time.time()

        if args.pre_pdb:
            structure_dir = args.pre_pdb
        else:
            logger.error("overall mode requires --pre_pdb")
            return

        # 2. feature calculation — only if PDB files exist
        pdb_list = [
            os.path.join(structure_dir, fn)
            for fn in os.listdir(structure_dir)
            if fn.lower().endswith(".pdb")
        ]

        if len(pdb_list) == 0:
            logger.error(f"No PDB structures found in: {structure_dir}")
            logger.error("Please provide --pre_pdb or run --mode fetch to fetch pdbs")
            return

        # 2. feature calculation
        return_mode = "frequency"
        feature_dict = calc_feature(pdb_list, return_mode)

        # 3. save overall
        overall_file = os.path.join(args.work_dir, "overall_feature.csv")
        overall.save_feature(overall_file, feature_dict)

        logger.info(f"Overall analysis finished in {int(time.time() - start)} seconds.")
        return

    # ------------------------------------------------------------------------------
    # MODE — LOCAL: Align local structure and analyze the features
    # ------------------------------------------------------------------------------
    elif args.mode == "local":
        start = time.time()
        if args.pre_pdb:
            structure_dir = args.pre_pdb
        else:
            logger.error("overall mode requires --pre_pdb")
            return

        if not (args.template_name and args.positions and args.dist1 and args.dist2):
            parser.error("Local mode requires --template_name --template_active_res --dist1 --dist2")

        pdb_list = [
            os.path.join(structure_dir, fn)
            for fn in os.listdir(structure_dir)
            if fn.lower().endswith(".pdb")
        ]

        if len(pdb_list) == 0:
            logger.error(f"No PDB structures found in: {structure_dir}")
            logger.error("Please provide --pre_pdb or run --mode fetch to fetch pdbs")
            return

        # 1. hydrophobic cluster calculation
        hydrophobic_feature = calc_local_hydrophobic(pdb_list)

        # 2. alignment
        align_file = align_structure(structure_dir)

        # 3. local output
        local_hpd_file = os.path.join(args.work_dir, "local_hydrophobic_feature.csv")
        local_aa_file = os.path.join(args.work_dir, "local_aa_feature.csv")

        # 5. run local analysis
        local.run(
            args.template_name,
            align_file,
            args.positions.split(","),
            hydrophobic_feature,
            pdb_list,
            local_hpd_file,
            local_aa_file,
            int(args.dist1),
            int(args.dist2)
        )

        logger.info(f"Local analysis finished in {int(time.time() - start)} seconds.")
        return

    # ------------------------------------------------------------------------------
    # MODE — VISUAL: Visualization of the four interactions in all pdbs
    # ------------------------------------------------------------------------------
    elif args.mode == "visual":
        if args.pre_pdb:
            structure_dir = args.pre_pdb
        else:
            logger.error("visual mode requires --pre_pdb")
            return

        internal.run(work_dir=args.work_dir, pdb_dir=structure_dir, pml=args.pml)

    # ------------------------------------------------------------------------------
    # MODE — surface: Surface charge analysis
    # ------------------------------------------------------------------------------
    elif args.mode == "surface":
        if args.pre_pdb:
            structure_dir = args.pre_pdb
        else:
            logger.error("visual mode requires --pre_pdb")
            return

        if not args.config_file:
            logger.error("visual mode requires --config_file")
            return

        surface.run(work_dir=args.work_dir, structure_dir=structure_dir, config_file=args.config_file)

    # ------------------------------------------------------------------------------
    # MODE — landscape: function landscape visualization
    # ------------------------------------------------------------------------------
    elif args.mode == "landscape":
        if args.pre_pdb:
            structure_dir = args.pre_pdb
        else:
            logger.error("visual mode requires --pre_pdb")
            return

        if not args.label_file:
            logger.error("visual mode requires --label_file")
            return

        if not args.template_name:
            logger.error("visual mode requires --template_name")
            return

        if not args.positions:
            logger.error("visual mode requires --positions")
            return

        if not args.label_pos_threshold:
            logger.error("visual mode requires --label_pos_threshold")
            return

        if not args.config_file:
            logger.error("visual mode requires --config_file")
            return

        landscape.run(
            work_dir=args.work_dir,
            pdb_dir=structure_dir,
            label_file=args.label_file,
            template_name=args.template_name,
            input_pos_list=args.positions,
            label_pos_threshold=args.label_pos_threshold,
            config_file=args.config_file,
        )

if __name__ == "__main__":
    main()
