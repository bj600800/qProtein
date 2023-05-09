"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/04/14

# Description: Structure quantitative characterization
# ------------------------------------------------------------------------------
"""
import os
import io
import csv

from tqdm import tqdm
import threading
import multiprocessing
from queue import Queue

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.DSSP import DSSP
import freesasa

from qprotein.database.sqlite3_builder import SqlBuilder
from qprotein.database.sqlite3_searcher import SqlSearch
from qprotein.utilities import logger
logger = logger.setup_log(name=__name__)


def get_dssp_cif(struct_path):
    parser = MMCIFParser()
    struct_file_name = os.path.split(struct_path)[1].split('.')[0]
    model = parser.get_structure(struct_file_name, struct_path)[0]
    # dssp = DSSP(model, struct_path, dssp="mkdssp")
    dssp = DSSP(model, struct_path, dssp=r"D:\subject\active\1-qProtein\tools\dssp\dssp-3.0.0-win32.exe")
    dssp_dat = [residue_dat for residue_dat in dssp]
    return dssp_dat


def get_dssp_pdb(struct_path):
    p = PDBParser()
    struct_file_name = os.path.split(struct_path)[1].split('.')[0]
    structure = p.get_structure(struct_file_name, struct_path)
    model = structure[0]
    dssp = DSSP(model, struct_path, dssp="mkdssp")
    dssp_dat = [residue_dat for residue_dat in dssp]
    return dssp_dat


def get_second_struct(dssp_dat):
    second_struct = [residue_dat[2].replace('-', 'C') for residue_dat in dssp_dat]
    ss = ''.join(second_struct).replace('G', 'H').replace('I', 'H')
    return ss


def get_ss_content(ss):
    ss_content = {}
    total_length = len(ss)
    ss_content['helix'] = ss.count('H')/total_length
    ss_content['sheet'] = ss.count('E')/total_length
    ss_content['loop'] = ss.count('C') / total_length
    ss_content['turn'] = ss.count('T')/total_length
    ss_content['bend'] = ss.count('S')/total_length
    ss_content['bridge'] = ss.count('B')/total_length
    return ss_content


def get_backbone_geometry(ss):
    ss_sequence = []
    i = 0
    while i < len(ss):
        c = ss[i]
        count = 1
        while i + count < len(ss) and ss[i + count] == c:
            count += 1
        ss_sequence.append(c + c * (count - 1))
        i += count


def get_Hbond(dssp_dat):
    """
    surface: Rsass > 0.05
    """
    aa_num = len(dssp_dat)
    hbond_num = 0
    hbond_energy = 0
    hbond_cols = [6, 8, 10, 12]
    hbond_energy_cols = [7, 9, 11, 13]
    for residue_dat in dssp_dat:
        for i, col in enumerate(hbond_cols):
            if residue_dat[col] != 0:
                hbond_num += 1
                hbond_energy += residue_dat[hbond_energy_cols[i]]
    hbond = {'hbond_density': hbond_num/aa_num, 'hbond_avg_energy': hbond_energy/aa_num}
    return hbond


def analyze_surface(dssp_dat):
    """
    surface: Rsass > 0.05
    the number of polar residues, non-polar residues, positive charge residues, negative charge residues.
    """
    # Define residue groups
    residue_groups = {
        'polar': ('S', 'T', 'Y', 'N', 'Q', 'C', 'K', 'R', 'H', 'D', 'E'),
        'apolar': ('G', 'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P'),
        'positive': ('K', 'R', 'H'),
        'negative': ('D', 'E')
        }
    # Filter surface residues
    surface_residues = [(dat[0], dat[1]) for dat in dssp_dat if dat[3] > 0.05]
    total_surface_residue = len(surface_residues)

    # Calculate proportions of different residue groups
    surface = {}
    for group, residues in residue_groups.items():
        count_group = len([res_site for res_site, res_type in surface_residues if res_type in residues])
        surface[group] = count_group / total_surface_residue
    return surface


def cif2pdb(struct_path):
    struct_name = os.path.splitext(os.path.basename(struct_path))[0]
    parser = MMCIFParser()
    structure = parser.get_structure(struct_name, struct_path)
    pdb_output = io.StringIO()
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(pdb_output)
    pdb_file = io.StringIO(pdb_output.getvalue())
    return pdb_file


def get_sasa(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure("pdb", pdb_file)
    result, sasa = freesasa.calcBioPDB(structure)
    # sasa-> {'Polar': 4657.344621590802, 'Apolar': 5662.672941379591}
    return sasa


def get_electrostatics():
    pass


def process_structure(struct_path):
    """
    Iterate every structure for process.
    """

    if struct_path.split('.')[-1] == 'cif':
        dssp_dat = get_dssp_cif(struct_path)
        pdb_file = cif2pdb(struct_path)
        sasa = get_sasa(pdb_file)
    else:
        dssp_dat = get_dssp_pdb(struct_path)
        sasa = get_sasa(struct_path)
    ss_content = get_ss_content(get_second_struct(dssp_dat))
    hbond = get_Hbond(dssp_dat)
    surface = analyze_surface(dssp_dat)
    return {'ss_content': ss_content, 'hbond': hbond, 'surface': surface, 'sasa': sasa}


class Producer(threading.Thread):
    def __init__(self, path_queue, results_queue, pbar):
        super(Producer, self).__init__()
        self.path_queue = path_queue
        self.results_queue = results_queue
        self.pbar = pbar

    def run(self):
        while True:
            if self.path_queue.empty():
                break
            struct_path = self.path_queue.get()
            struct_name = os.path.split(struct_path)[1].split('.')[0]
            results = process_structure(struct_path)
            self.results_queue.put((struct_name, results))
            if self.pbar:
                self.pbar.update(1)



class Consumer(threading.Thread):
    def __init__(self, output_file, path_queue, results_queue, pbar):
        super(Consumer, self).__init__()
        self.path_queue = path_queue
        self.results_queue = results_queue
        self.output_file = output_file
        self.pbar = pbar

    def insert_sql(self):
        # Todo
        pass

    def run(self):
        with open(self.output_file, 'a+', newline='') as csvfile:
            while True:
                if self.results_queue.empty():
                    break
                (struct_name, results) = self.results_queue.get()
                writer = csv.writer(csvfile)
                writer.writerow([struct_name] + [str(value) for subdict in results.values() for value in subdict.values()])
                # insert_many()
                if self.pbar:
                    self.pbar.update(1)


def get_process_num():
    core_num = multiprocessing.cpu_count()
    return core_num


def analyze(struct_dir, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['', 'helix', 'sheet', 'loop', 'turn', 'bend', 'bridge', 'hbond_density',
                         'hbond_avg_energy', 'apolar', 'polar', 'positive', 'negative',
                         'polar_area', 'apolar_area'])

    path_queue = Queue(maxsize=0)
    results_queue = Queue(maxsize=0)
    producer_thread = get_process_num()
    logger.info(f"Start {producer_thread} threads for analyzing the structures")

    struct_path_list = [os.path.join(struct_dir, i) for i in os.listdir(struct_dir) if os.path.splitext(i)[1] == '.cif'
                        and os.path.getsize(os.path.join(struct_dir, i)) > 0]
    if struct_path_list:
        for path in struct_path_list:
            path_queue.put(path)
    else:
        logger.info("No structure needs to analyze")
        return

    producer_thread_list = []
    with tqdm(total=len(struct_path_list), desc="Structure analyzing progress") as pbar:
        for i in range(producer_thread):
            t = Producer(path_queue=path_queue, results_queue=results_queue, pbar=pbar)
            producer_thread_list.append(t)
            t.start()
        for i in producer_thread_list:
            i.join()
    logger.info('Start to save the computing results')

    with tqdm(total=results_queue.qsize(), desc="Results saving progress") as pbar:
        consumer = Consumer(path_queue=path_queue, results_queue=results_queue,
                            output_file=output_file, pbar=pbar)
        consumer.start()
        consumer.join()
    logger.info('Structure analyzing mission finished!')

