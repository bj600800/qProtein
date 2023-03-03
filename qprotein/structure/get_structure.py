"""
Get Structure from 3D beacons
Prioritize protein structure quality:
PDB -> AlphaFold
"""

import threading  # 导入threading模块
from queue import Queue  # 导入queue模块
import time  # 导入time模块
from requests.adapters import HTTPAdapter, Retry
import requests
import json
import os


class GetStruct:
    def __init__(self):
        self.base_url = 'https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/v1/uniprot/summary/'
        retries = Retry(total=5, backoff_factor=0.25,
                        status_forcelist=[429, 500, 502, 503, 504])
        self.session = requests.Session()
        self.session.mount("https://", HTTPAdapter(max_retries=retries))

    def check_response(self, response):
        try:
            response.raise_for_status()
        except requests.HTTPError:
            return response.status_code

    def get_cif_url(self, uniID):
        """
        根据uniprotID爬取蛋白质结构链接
        """
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.25 Safari/537.36 Core/1.70.3861.400 QQBrowser/10.7.4313.400"
        }
        base_url = self.base_url

        para_struct = {}
        time1 = time.time()
        while True:
            # request for json
            try:
                r = self.session.get(url=base_url+uniID +
                                     '.json', headers=headers)

                res_json = json.loads(r.content)
                if res_json['uniprot_entry']['sequence_length']:
                    sequence_length = res_json['uniprot_entry']['sequence_length']
                    para_struct['sequence_length'] = sequence_length

                for structure in res_json['structures']:
                    provider = structure['provider']
                    para_struct['provider'] = provider

                    # parser for structure URL type
                    if provider == 'PDBe':
                        para_struct['identifier'] = structure['model_identifier']
                        para_struct['resolution'] = structure['resolution']
                        para_struct['cif_url'] = structure['model_url']
                        time2 = time.time()-time1
                        print("Get {uniID}(resolution: {resolution}, sequence_length: {sequence_length}) url {cif_url} from {provider} in {time2} seconds".format(
                            resolution=para_struct['resolution'], uniID=uniID, cif_url=para_struct['cif_url'], provider=para_struct['provider'], sequence_length=para_struct['sequence_length'], time2=time2))
                        return para_struct['cif_url']

                    elif provider == 'AlphaFold DB':
                        para_struct['identifier'] = structure['model_identifier']
                        para_struct['cif_url'] = structure['model_url']
                        time2 = time.time()-time1
                        print("Get {uniID}(sequence_length: {sequence_length}) url {cif_url} from {provider} in {time2} seconds".format(
                            uniID=uniID, cif_url=para_struct['cif_url'], provider=para_struct['provider'], sequence_length=para_struct['sequence_length'], time2=time2))
                        return para_struct['cif_url']
            except:
                print('error during requesting structrue URL: ', uniID)
                print('-'*30)
                break

    def get_structure(self, url):
        """
        下载蛋白质结构
        """
        if not os.path.exists('str/'):
            os.mkdir('str/')

        time1 = time.time()
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.25 Safari/537.36 Core/1.70.3861.400 QQBrowser/10.7.4313.400"
        }
        while True:
            try:
                cif = self.session.get(url, headers=headers)
                cif_str = cif.content.decode('utf-8')

                structure_name = url.split('/')[-1][:-4]
                file_path = os.path.join('str', structure_name)

                with open(file_path+'.cif', 'w') as f:
                    f.writelines(cif_str)
                print('Structure downloaded successfully!')
                break

            except:
                print('error during downloading structure url: ', url)
                break

        time2 = time.time()-time1
        print('{url} download time: {time2}'.format(url=url, time2=time2))
        print('-'*30)

    def main(self, uniID):
        start_time = time.time()

        if not os.path.exists('E:\pymol\Lib\site-packages\pmg_tk\startup\pyMulstruct\server\src\strucutre_check.txt'):
            open('E:\pymol\Lib\site-packages\pmg_tk\startup\pyMulstruct\server\src\strucutre_check.txt', 'w')

        with open('E:\pymol\Lib\site-packages\pmg_tk\startup\pyMulstruct\server\src\list.txt', 'r') as f1:
            meta_list = [i.strip() for i in f1.readlines()]
            print(meta_list)
            with open('E:\pymol\Lib\site-packages\pmg_tk\startup\pyMulstruct\server\src\strucutre_check.txt', 'r') as f2:
                log_list = [i.strip() for i in f2.readlines()]
                print(log_list)
                continue_list = [i for i in meta_list if i not in log_list]

            for uniID in continue_list:
                url = self.get_cif_url(uniID)
                if url:
                    self.get_structure(url)
                with open('E:\pymol\Lib\site-packages\pmg_tk\startup\pyMulstruct\server\src\strucutre_check.txt', 'a') as wr:
                    wr.write(uniID+'\n')

        print("last time: {} s".format(time.time()-start_time))


if __name__ == "__main__":
    run = GetStruct()
    run.main()
