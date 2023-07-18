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
from Bio import BiopythonWarning
import warnings
warnings.simplefilter("ignore", BiopythonWarning)

from qprotein.utilities import logger
logger = logger.setup_log(name=__name__)

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

    struct_path_list = [os.path.join(struct_dir, i) for i in os.listdir(struct_dir)]
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

if __name__ == '__main__':
    struct_dir = r"D:\subject\active\1-qProtein\data\enzymes\cazy_endo-1_4-beta-xylanase\classification\thermophilic"
    output_file = r"D:\subject\active\1-qProtein\data\enzymes\cazy_endo-1_4-beta-xylanase\classification\thermophilic.out"
    analyze(struct_dir=struct_dir, output_file=output_file)
