"""
Get Structure from 3D beacons
Prioritize protein structure quality:
PDB -> AlphaFold
"""

import json
import os
import sqlite3
import threading
from queue import Queue

import requests


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise


def sql_parse(sql_db):
    connect = sqlite3.connect(sql_db)
    cursor = connect.cursor()
    sql = "SELECT query_name, cazy_family, sprot_acc, trembl_acc FROM summary WHERE trembl100_acc !='' or sprot100_acc !=''"
    cursor.execute("begin")
    cursor.execute(sql)
    sql_rows = cursor.fetchall()
    print('Get targets: ', len(sql_rows))

    dict_uniprot_ids = {}
    for i in sql_rows:
        if i[1] != None:
            dict_uniprot_ids[i[0] + '_sprot_' + i[1]] = i[1]
        elif i[2] != None:
            dict_uniprot_ids[i[0] + '_trembl_' + i[2]] = i[2]
    print('SQL search complete!')
    return dict_uniprot_ids


class Producer(threading.Thread):
    def __init__(self, uniID_queue, struct_queue):
        super(Producer, self).__init__()
        self.POLLING_INTERVAL = 1
        self.base_url = 'https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/'
        self.uniID_queue = uniID_queue
        self.struct_queue = struct_queue

    def run(self) -> None:
        while True:
            if self.uniID_queue.empty():
                break
            file_name, uniID = self.uniID_queue.get()
            self.get_cif_url(file_name, uniID)

    def get_cif_url(self, file_name, uniID):
        """
        根据uniprotID爬取蛋白质结构链接
        """
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.25 Safari/537.36 Core/1.70.3861.400 QQBrowser/10.7.4313.400"
            }
        base_url = self.base_url
        print('Get uniID: ', uniID)

        para_struct = {}
        try:
            r = requests.get(url=base_url + uniID + '.json', headers=headers)
            check_response(r)
            ret_json = json.loads(r.content)
            print(ret_json)

            if ret_json['uniprot_entry']['sequence_length']:
                sequence_length = ret_json['uniprot_entry']['sequence_length']
                para_struct['sequence_length'] = sequence_length

            for structure in ret_json['structures']:
                provider = structure['summary']['provider']
                para_struct['provider'] = provider

                if provider == 'PDBe':

                    if structure['summary']['model_identifier']:
                        identifier = structure['summary']['model_identifier']
                        para_struct['identifier'] = identifier
                    else:
                        para_struct['identifier'] = ''

                    if structure['resolution']:
                        resolution = structure['summary']['resolution']
                        para_struct['resolution'] = resolution
                    else:
                        para_struct['resolution'] = ''

                    if structure['model_url']:
                        cif_url = structure['summary']['model_url']
                        print(cif_url)
                        print(cif_url)
                        para_struct['cif_url'] = cif_url
                        self.struct_queue.put((file_name, cif_url))
                        print(file_name, cif_url)
                    else:
                        para_struct['cif_url'] = ''

                    # Queue队列的put方法用于向Queue队列中放置元素，由于Queue是先进先出队列，所以先被Put的URL也就会被先get出来。

                    # print(
                    #     "Thread {thread} get {identifier}(resolution: {resolution}, sequence_length: {sequence_length}) url {cif_url} from {provider}".format(
                    #         thread=thread, resolution=para_struct['resolution'], identifier=para_struct['identifier'],
                    #         cif_url=para_struct['cif_url'], provider=para_struct['provider'],
                    #         sequence_length=para_struct['sequence_length']))
                    break

                elif provider == 'AlphaFold DB':

                    if structure['summary']['model_identifier']:
                        identifier = structure['summary']['model_identifier']
                        para_struct['identifier'] = identifier
                    else:
                        para_struct['identifier'] = ''

                    if structure['summary']['model_url']:
                        cif_url = structure['summary']['model_url']
                        print(cif_url)
                        para_struct['cif_url'] = cif_url
                        self.struct_queue.put((file_name, cif_url))
                        print(file_name, cif_url)
                    else:
                        para_struct['cif_url'] = ''

                    # Queue队列的put方法用于向Queue队列中放置元素，由于Queue是先进先出队列，所以先被Put的URL也就会被先get出来。

                    # print("get {identifier}(sequence_length: {sequence_length}) url {cif_url} from {provider}".format(
                    #     identifier=para_struct['identifier'], cif_url=para_struct['cif_url'],
                    #     provider=para_struct['provider'], sequence_length=para_struct['sequence_length']))
                    break
        except:
            print('url error')
            pass


class Consumer(threading.Thread):
    def __init__(self, uniID_queue, struct_queue):
        super(Consumer, self).__init__()
        self.POLLING_INTERVAL = 1
        self.uniID_queue = uniID_queue
        self.struct_queue = struct_queue

    def run(self):
        while True:
            if self.uniID_queue.empty() and self.struct_queue.empty():
                print('empty')
                break
            file_name, url = self.struct_queue.get()
            print(file_name, url)
            self.get_struct(file_name, url)

    def get_struct(self, file_name, url):
        """
        下载蛋白质结构
        """
        print("Downloaded {url}".format(url=url))  # 打印线程thread和被爬取结构的url
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.71 Safari/537.36"
            }
        try:
            cif = requests.get(url, headers=headers)
            check_response(cif)
            cif_str = cif.content.decode('utf-8')
            print('Downloading ', url)
            dir_path = r'D:\subject\active\1-PyMulstruct\data\tibet\str1'
            if not os.path.exists(dir_path):
                os.mkdir(dir_path)
            file_path = os.path.join(dir_path, file_name)
            with open(file_path + '.cif', 'w') as f:
                f.writelines(cif_str)


        except requests.HTTPError:
            print('Error: ', file_name)
            raise


def main(sql_db):
    # 用Queue构造一个大小为1000的线程安全的先进先出队列
    uniID_queue = Queue(1000)
    struct_queue = Queue(1000)
    dict_uniprot_ids = sql_parse(sql_db)
    # put uniID_queue
    for file_name, uniID in dict_uniprot_ids.items():
        print("Put queue ", uniID)
        uniID_queue.put((file_name, uniID))

    # define 5 producer
    producer_thread = 5
    producer_thread_list = []
    for i in range(producer_thread):
        t = Producer(uniID_queue=uniID_queue, struct_queue=struct_queue)
        print(t)
        producer_thread_list.append(t)
        t.start()

    # define 3 consumer
    consumer_thread = 5
    consumer_thread_list = []
    for i in range(consumer_thread):
        t = Consumer(uniID_queue=uniID_queue, struct_queue=struct_queue)
        print(t)
        consumer_thread_list.append(t)
        t.daemon = True
        t.start()

    for t in producer_thread_list:
        t.join()


if __name__ == "__main__":
    sql_db = r'D:\subject\active\1-PyMulstruct\data\tibet\summary.db'
    main(sql_db)
