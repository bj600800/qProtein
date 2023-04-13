import json
import os
import threading
from queue import Queue

import requests
import multiprocessing
from qprotein.database.sqlite3_searcher import SqlSearch
from qprotein.structure.StructureSegmenter import Segmenter
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


class TargetNameRetriever(SqlSearch):
    """
    Get target candidates from sql database
    """

    def __init__(self, summary_results_db):
        super(TargetNameRetriever, self).__init__(sql_db=summary_results_db)

    def get_sql_results(self):
        cursor = self.connect_sql(sql_db=sql_db)
        sql_cmd = "SELECT query_name, seq_length, sprot_acc, sprot_ident, sprot_cover, sprot_match_length," \
                  "sprot_query_start, sprot_query_end, sprot_subject_start, sprot_subject_end," \
                  "trembl_acc, trembl_ident, trembl_cover, trembl_match_length, trembl_query_start," \
                  "trembl_query_end, trembl_subject_start, trembl_subject_end  " \
                  "FROM results_summary WHERE sprot_acc != '' or trembl_acc != ''"
        results_list = self.fetch_results(cursor=cursor, sql_cmd=sql_cmd)
        return results_list

    def filter_results(self):
        # concerned features: identity, coverage
        results_list = self.get_sql_results()
        results_filtered = []
        for result in results_list:
            if all(bool(i) for i in (result[2], result[10])):
                sprot_quality = (result[3], result[4])
                trembl_quality = (result[11], result[12])
                if sum(sprot_quality) >= sum(trembl_quality):
                    results_filtered.append(result[:10])
                else:
                    results_filtered.append(result[:2] + result[10:])
            elif result[2] is not None:
                results_filtered.append(result[:10])
            elif result[10] is not None:
                results_filtered.append(result[:2] + result[10:])
        return results_filtered


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError as e:
        logger.debug(e)


class Producer(threading.Thread):
    """
    Produce url of structures from AFDB
    """

    def __init__(self, query_info_queue, structure_info_queue):
        super(Producer, self).__init__()
        self.base_url = "https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/"
        self.query_info_queue = query_info_queue
        self.structure_info_queue = structure_info_queue

    def get_cif_url(self, uniprot_acc):
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.25 Safari/537.36 Core/1.70.3861.400 QQBrowser/10.7.4313.400"
            }
        base_url = self.base_url
        try:
            r = requests.get(url=base_url + uniprot_acc + '.json?provider=alphafold', headers=headers)
            check_response(r)
            ret_json = json.loads(r.content)
            return ret_json

        except requests.HTTPError as e:
            logger.debug(uniprot_acc + ": failed for " + e)

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
            # protein_info = ('Sum14_k97_629838_2_1', 383, 'A0A1M3MXV5', 92.4, 100.0, 383, 1, 383, 77, 459)


class Consumer(threading.Thread):
    """
    get structures and cut them off
    """

    def __init__(self, query_info_queue, structure_info_queue, dir_path):
        super(Consumer, self).__init__()
        self.query_info_queue = query_info_queue
        self.structure_info_queue = structure_info_queue
        self.dir_path = dir_path

    @staticmethod
    def request_structure(uniprot_acc, url):
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.71 Safari/537.36",
            "From": "bj600800@gmail.com"  # ALLWAYS TELL WHO YOU ARE
            }
        try:
            cif = requests.get(url, headers=headers)
            check_response(cif)
            cif_str = cif.content.decode('utf-8')
            return cif_str

        except requests.HTTPError as e:
            logger.debug(uniprot_acc + ": failed for " + e)
            raise

    def map_structure(self, cif_str, query_info):
        query_name, query_length, uniprot_acc, identity, coverage, query_match_length, query_start_residue, \
        query_end_residue, subject_start_residue, subject_end_residue = query_info

        write_cif_path_kw = [self.dir_path, query_name, uniprot_acc]

        segmenter = Segmenter(cif_string=cif_str, query_name=query_name, subject_start_residue=subject_start_residue,
                              subject_end_residue=subject_end_residue, query_start_residue=query_start_residue,
                              query_end_residue=query_end_residue, query_match_length=query_match_length,
                              query_length=query_length, write_cif_path_kw=write_cif_path_kw
                              )
        write_cif_path = segmenter.run()
        return write_cif_path

    def insert_sql(self, info):
        pass

    def run(self):
        while True:
            if self.structure_info_queue.empty() and self.query_info_queue.empty():
                logger.info('Structure mapping mission complete!')
                break
            query_info, structure_info = self.structure_info_queue.get()
            uniprot_acc = structure_info["uniprot_acc"]
            url = structure_info["url"]
            cif_str = self.request_structure(uniprot_acc=uniprot_acc, url=url)
            self.map_structure(cif_str=cif_str, query_info=query_info)


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


def main(summary_results_db, dir_path):
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)

    query_info_queue = Queue(maxsize=0)
    structure_info_queue = Queue(maxsize=0)
    target_handler = TargetNameRetriever(summary_results_db=summary_results_db)
    query_info = target_handler.filter_results()
    logger.info(f'Found {len(query_info)} targets to map structure.')
    for info in query_info:
        query_info_queue.put(info)

    producer_thread, consumer_thread = get_producer_comsumer_num()
    logger.info(f'We have {producer_thread} thread(s) for producer and {consumer_thread} thread(s) for consumer.')

    producer_thread_list = []
    for i in range(producer_thread):
        t = Producer(query_info_queue=query_info_queue, structure_info_queue=structure_info_queue)
        producer_thread_list.append(t)
        t.start()

    consumer_thread_list = []
    for i in range(consumer_thread):
        t = Consumer(query_info_queue=query_info_queue, structure_info_queue=structure_info_queue, dir_path=dir_path)
        consumer_thread_list.append(t)
        t.daemon = True
        t.start()

    for t in consumer_thread_list:
        t.join()

    for t in producer_thread_list:
        t.join()


if __name__ == "__main__":
    sql_db = r'D:\subject\active\1-qProtein\data\manure\qprotein_results.db'
    dir_path = r'D:\subject\active\1-qProtein\data\manure\structure'
    main(sql_db, dir_path)