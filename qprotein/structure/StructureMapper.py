"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/04/14

# Description: Mapping Structure for the target sequence
# ------------------------------------------------------------------------------
"""
import json
import os
import threading
import multiprocessing
from queue import Queue
import time
import random
from tqdm import tqdm
import requests
from requests.adapters import Retry
from pdbecif.mmcif_io import CifFileWriter

from qprotein.database.sqlite3_searcher import SqlSearch
from qprotein.structure.StructureSegmenter import Segmenter
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


class TargetNameRetriever(SqlSearch):
    """
    Get target candidates from sql database
    """

    def __init__(self, summary_results_db, structure_dir_path, log_file):
        super(TargetNameRetriever, self).__init__()
        self.sql_db = summary_results_db
        self.structure_dir_path = structure_dir_path
        self.log_file = log_file
    
    def get_sql_results(self):
        cursor = self.connect_sql(sql_db=self.sql_db)
        sql_cmd = "SELECT query_name, seq_length, sprot_acc, sprot_subject_start, sprot_subject_end," \
                  "trembl_acc, trembl_subject_start, trembl_subject_end  " \
                  "FROM results_summary WHERE sprot_acc != '' or trembl_acc != ''"

        results_list = self.fetch_results(cursor=cursor, sql_cmd=sql_cmd)
        logger.info(f"Get the total number of targets from SQL: {len(results_list)}")
        return results_list

    def get_exist_query(self):
        exist_structure = [file for path in self.structure_dir_path for file in os.listdir(path)
                           if os.path.isfile(os.path.join(path, file)) and os.path.getsize(os.path.join(path, file)) > 0]
        exist_query_name = [item.split('#')[1] for item in exist_structure]
        return exist_query_name

    def get_NO_structure_query(self):
        if os.path.exists(self.log_file):
            with open(self.log_file, 'r') as rf:
                invalid_query_name = [i.rstrip() for i in rf.readlines()]
            return invalid_query_name
        else:
            return []

    def filter_results(self):
        exist_query_name = self.get_exist_query()
        logger.info(f'Found {len(exist_query_name)} structure(s)')
        invalid_query_name = self.get_NO_structure_query()
        logger.info(f'Found {len(invalid_query_name)} non-structural targets')

        results_list = [result for result in tqdm(self.get_sql_results(), desc='Filter query targets', unit='step')
                        if result[0] not in exist_query_name
                        if result[2] not in invalid_query_name
                        and result[10] not in invalid_query_name]
        query_info = []
        for result in results_list:
            if all(bool(i) for i in (result[2], result[10])):
                sprot_quality = (result[3], result[4])
                trembl_quality = (result[11], result[12])
                if sum(sprot_quality) >= sum(trembl_quality):
                    query_info.append(result[:10] + ("sprot",))
                else:
                    query_info.append(result[:2] + result[10:] + ("trembl",))
            elif result[2] is not None:
                query_info.append(result[:10] + ("sprot",))
            elif result[10] is not None:
                query_info.append(result[:2] + result[10:] + ("trembl",))
        return query_info


def start_request_session():
    retries = Retry(total=3, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504])
    session = requests.Session()
    session.mount('http://', requests.adapters.HTTPAdapter(max_retries=retries))
    session.mount('https://', requests.adapters.HTTPAdapter(max_retries=retries))
    return session


class Producer(threading.Thread):
    """
    Produce url of structures from AFDB
    """

    def __init__(self, query_info_queue, structure_info_queue, log_file):
        super(Producer, self).__init__()
        self.base_url = "https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/"
        self.query_info_queue = query_info_queue
        self.structure_info_queue = structure_info_queue
        self.log_file = log_file
        self.exist_invalid_name = []

    def get_existed_invalid_uniprot(self):
        if os.path.exists(self.log_file):
            with open(self.log_file, 'r') as rf:
                self.exist_invalid_name = [i.strip() for i in rf.readlines()]

    def check_response(self, response, uniprot_acc, log_file):
        try:
            response.raise_for_status()
        except requests.HTTPError:
            if response.status_code == 404:
                self.get_existed_invalid_uniprot()
                if uniprot_acc not in self.exist_invalid_name:
                    with open(log_file, 'a+') as wf:
                        wf.write(uniprot_acc+'\n')
            logger.debug(f"Not found structure for {uniprot_acc}")

    def get_cif_url(self, uniprot_acc):
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.25 Safari/537.36 Core/1.70.3861.400 QQBrowser/10.7.4313.400",
            "From": "bj600800@gmail.com"  # ALLWAYS TELL WHO YOU ARE
            }
        base_url = self.base_url
        session = start_request_session()
        r = session.get(url=base_url + uniprot_acc + '.json?provider=alphafold', headers=headers)
        time.sleep(random.uniform(1, 3))  # DO NOT ABUSE THE SERVER
        self.check_response(r, uniprot_acc=uniprot_acc, log_file=self.log_file)
        ret_json = json.loads(r.content)
        return ret_json

    @staticmethod
    def parse_json(ret_json, uniprot_acc):
        structure_info = {"uniprot_acc": uniprot_acc}
        uniprot_entry = ret_json["uniprot_entry"]
        if uniprot_entry["ac"] == uniprot_acc:
            structure_info["sequence_length"] = uniprot_entry.get("sequence_length", "")
            for structure in ret_json["structures"]:
                summary = structure["summary"]
                structure_info["plddt"] = summary.get("confidence_avg_local_score", "")
                structure_info["url"] = summary.get("model_url", "")
            return structure_info

    def put_structure_info(self, query_info):
        uniprot_acc = query_info[2]
        ret_json = self.get_cif_url(uniprot_acc=uniprot_acc)
        if ret_json:
            structure_info = self.parse_json(ret_json=ret_json, uniprot_acc=uniprot_acc)
            if structure_info["url"]:
                all_info = (query_info, structure_info)
                self.structure_info_queue.put(all_info)

    def run(self):
        while True:
            if self.query_info_queue.empty():
                break
            query_info = self.query_info_queue.get()
            self.put_structure_info(query_info)


class Consumer(threading.Thread):
    """
    get structures and cut them off
    """

    def __init__(self, query_info_queue, structure_info_queue, structure_dir_path):
        super(Consumer, self).__init__()
        self.query_info_queue = query_info_queue
        self.structure_info_queue = structure_info_queue
        self.structure_dir_path = structure_dir_path

    @staticmethod
    def check_response(response):
        try:
            response.raise_for_status()
        except requests.HTTPError as e:
            logger.debug(e)

    def request_structure(self, url):
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.71 Safari/537.36",
            "From": "bj600800@gmail.com"  # ALLWAYS TELL the SERVER WHO YOU ARE
            }
        session = start_request_session()
        response = session.get(url, headers=headers)
        self.check_response(response=response)
        time.sleep(random.uniform(1, 3))  # DO NOT ABUSE THE SERVER
        cif_str = response.content.decode('utf-8')
        return cif_str

    def write_cif(self, mean_plddt, output_cif, uniprot_class, query_name, uniprot_acc):
        high_plddt_dir_path, confident_plddt_dir_path = self.structure_dir_path

        if mean_plddt >= 90:
            write_cif_path = os.path.join(high_plddt_dir_path,
                                          '#'.join([uniprot_class, query_name, uniprot_acc, str(mean_plddt)])
                                          ) + ".cif"
            writer = CifFileWriter(write_cif_path)
            writer.write(output_cif)
            return write_cif_path
        elif 70 <= mean_plddt < 90:
            write_cif_path = os.path.join(confident_plddt_dir_path,
                                          '#'.join([uniprot_class, query_name, uniprot_acc, str(mean_plddt)])
                                          ) + ".cif"
            writer = CifFileWriter(write_cif_path)
            writer.write(output_cif)
            return write_cif_path
        else:
            return f"Skipped for low average pLDDT -> {mean_plddt}"

    def map_structure(self, cif_str, query_info):
        query_name, query_length, uniprot_acc, identity, coverage, query_match_length, query_start_residue, \
            query_end_residue, subject_start_residue, subject_end_residue, uniprot_class = query_info

        write_cif_path_kw = [self.structure_dir_path, query_name, uniprot_acc]
        segmenter = Segmenter(cif_string=cif_str, query_name=query_name, subject_start_residue=subject_start_residue,
                              subject_end_residue=subject_end_residue, query_start_residue=query_start_residue,
                              query_end_residue=query_end_residue, query_match_length=query_match_length,
                              query_length=query_length, write_cif_path_kw=write_cif_path_kw, uniprot_class=uniprot_class
                              )
        try:
            output_cif, mean_plddt = segmenter.run()
            write_cif_path = self.write_cif(output_cif=output_cif, mean_plddt=mean_plddt, uniprot_class=uniprot_class,
                                            query_name=query_name, uniprot_acc=uniprot_acc)
            return write_cif_path
        except:
            return

    def insert_sql(self, info):
        pass

    def run(self):
        while True:
            if self.structure_info_queue.empty() and self.query_info_queue.empty():
                logger.info('Structure mapping mission complete!')
                break
            query_info, structure_info = self.structure_info_queue.get()
            url = structure_info["url"]
            cif_str = self.request_structure(url=url)
            write_cif_path = self.map_structure(cif_str=cif_str, query_info=query_info)
            logger.info(f'Write {query_info[0]} structure: {write_cif_path}')



def get_cpu_core():
    core_num = multiprocessing.cpu_count()
    return core_num


def get_producer_comsumer_num():
    core_num = get_cpu_core()
    for producer_num in range(core_num):
        consumer_num = core_num - producer_num
        if producer_num >= 2 * consumer_num:
            return producer_num, consumer_num
    return 1, 1


def map_multiprocess(task_dir):
    summary_results_db = os.path.join(task_dir, "qprotein_results.db")
    log_file = os.path.join(task_dir, "404NotFoundURL.txt")
    structure_dir = ("high_structure", "confident_structure")
    structure_dir_path = [os.path.join(task_dir, i) for i in structure_dir]
    for path in structure_dir_path:
        if not os.path.exists(path):
            os.mkdir(path)
    query_info_queue = Queue(maxsize=0)
    structure_info_queue = Queue(maxsize=0)
    target_handler = TargetNameRetriever(summary_results_db=summary_results_db, structure_dir_path=structure_dir_path,
                                         log_file=log_file)
    query_info = target_handler.filter_results()
    logger.info(f'Found {len(query_info)} targets to map structure.')
    if query_info:
        for info in query_info:
            query_info_queue.put(info)
        producer_thread, consumer_thread = get_producer_comsumer_num()
        logger.info(f'We have {producer_thread} thread(s) for producer and {consumer_thread} thread(s) for consumer.')

        producer_thread_list = []
        for i in range(producer_thread):
            t = Producer(query_info_queue=query_info_queue, structure_info_queue=structure_info_queue, log_file=log_file)
            producer_thread_list.append(t)
            t.start()

        consumer_thread_list = []
        for i in range(consumer_thread):
            t = Consumer(query_info_queue=query_info_queue, structure_info_queue=structure_info_queue,
                         structure_dir_path=structure_dir_path)
            consumer_thread_list.append(t)
            t.daemon = True
            t.start()

        for t in consumer_thread_list:
            t.join()

        for t in producer_thread_list:
            t.join()

        total_structure = 0
        for dir_path in structure_dir_path:
            total_structure += len(os.listdir(dir_path))
        logger.info(f'qProtein collected {total_structure} structures in total')
    else:
        logger.info("No query sequence needs to map structure")


if __name__ == '__main__':
    task_dir = r"C:\Users\bj600\Desktop\test"
    map_multiprocess(task_dir)
