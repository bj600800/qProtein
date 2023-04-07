"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: obsoleted
  Copyright: Adapted from 'https://www.uniprot.org/help/id_mapping'
# ------------------------------------------------------------------------------
"""

import re
import time
import json
import zlib
import numpy as np
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry


class IDmapping:
    def __init__(self):
        self.POLLING_INTERVAL = 1
        self.API_URL = "https://rest.uniprot.org"
        retries = Retry(total=5, backoff_factor=0.25,
                        status_forcelist=[500, 502, 503, 504])
        self.session = requests.Session()
        self.session.mount("https://", HTTPAdapter(max_retries=retries))

    def check_response(self, response):
        try:
            response.raise_for_status()
        except requests.HTTPError:
            print(response.json()['messages'])
            raise

    def submit_id_mapping(self, input_id):
        from_db = ''
        to_db = ''

        if re.match('^[A-Z]{3}.+', input_id):
            from_db = 'EMBL-GenBank-DDBJ_CDS'  # 多个来源，判断refseq类型
            to_db = 'UniProtKB'
        elif re.match('^NP|^XP.+', input_id):
            from_db = 'RefSeq_Protein'
            to_db = 'UUniProtKB'

        elif re.match('^[A-Z]\d+', input_id):  # used for test
            from_db = 'UniProtKB_AC-ID'
            to_db = 'RefSeq_Protein'
        try:
            request = requests.post(
                f"{self.API_URL}/idmapping/run",
                data={"from": from_db, "to": to_db, "ids": input_id},
            )
            self.check_response(request)
            # print('submit_id_mapping succeed')
            return request.json()["jobId"]

        except:
            print('id {input_id} has an error'.format(input_id=input_id))

    def get_next_link(self, headers):
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def check_id_mapping_results_ready(self, job_id):
        while True:
            r = self.session.get(f"{self.API_URL}/idmapping/status/{job_id}")
            try:
                self.check_response(r)
                j = r.json()
                if "jobStatus" in j:
                    if j["jobStatus"] == "RUNNING":
                        print(f"Retrying in {self.POLLING_INTERVAL}s")
                        time.sleep(self.POLLING_INTERVAL)
                    else:
                        raise Exception(r["jobStatus"])
                else:
                    return bool(j["results"] or j["failedIds"])
            except:
                # print('check_id_mapping_results_ready error')
                break

    def get_batch(self, batch_response, file_format, compressed):
        batch_url = self.get_next_link(batch_response.headers)
        while batch_url:
            batch_response = self.session.get(batch_url)
            batch_response.raise_for_status()
            yield self.decode_results(batch_response, file_format, compressed)
            batch_url = self.get_next_link(batch_response.headers)

    def get_id_mapping_results_link(self, job_id):
        url = f"{self.API_URL}/idmapping/details/{job_id}"
        try:
            request = self.session.get(url)
            self.check_response(request)
            return request.json()["redirectURL"]
        except:
            pass

    def decode_results(self, response, file_format, compressed):
        if compressed:
            decompressed = zlib.decompress(
                response.content, 16 + zlib.MAX_WBITS)
            if file_format == "json":
                j = json.loads(decompressed.decode("utf-8"))
                return j
            elif file_format == "tsv":
                return [line for line in decompressed.decode("utf-8").split("\n") if line]
            elif file_format == "xlsx":
                return [decompressed]
            elif file_format == "xml":
                return [decompressed.decode("utf-8")]
            else:
                return decompressed.decode("utf-8")
        elif file_format == "json":
            return response.json()
        elif file_format == "tsv":
            return [line for line in response.text.split("\n") if line]
        elif file_format == "xlsx":
            return [response.content]
        elif file_format == "xml":
            return [response.text]
        return response.text

    def get_id_mapping_results_search(self, url):
        try:
            parsed = urlparse(url)
            query = parse_qs(parsed.query)
            file_format = query["format"][0] if "format" in query else "json"
            if "size" in query:
                size = int(query["size"][0])
            else:
                size = 500
                query["size"] = size
            compressed = (
                query["compressed"][0].lower(
                ) == "true" if "compressed" in query else False
            )
            parsed = parsed._replace(query=urlencode(query, doseq=True))
            url = parsed.geturl()
            request = self.session.get(url)
            self.check_response(request)
            results = self.decode_results(request, file_format, compressed)
            return results
        except:
            pass

    def run(self, input_id):
        job_id = self.submit_id_mapping(input_id=input_id)
        if self.check_id_mapping_results_ready(job_id):
            link = self.get_id_mapping_results_link(job_id)
            results = self.get_id_mapping_results_search(link)
            print(results['results'][0]['to'])
            return results  # 重新构造一个字典，只保留需要的信息，如序列，功能域注释等。
            # results['results'][0]['to']['sequence']['value']

        # Equivalently using the stream endpoint which is more demanding
        # on the API and so is less stable:
        # results = get_id_mapping_results_stream(link)

        # {'results': [{'from': 'P05067', 'to': 'CHEMBL2487'}], 'failedIds': ['P12345']}


if __name__ == '__main__':
    run = IDmapping()
    id_list = ['ASI30444.1']
    info = []
    for id in id_list:
        info.append(run.run(input_id=id))
        if len(info) % 2 == 0:
            np.save('../../tools/dict.npy', info)
            print('write {times}'.format(times=len(info)))
    if info[0]:
        print(info[0]['results'])
    # dict_keys(['entryType', 'primaryAccession', 'uniProtkbId', 'entryAudit', 'annotationScore', 'organism',
    # 'proteinExistence', 'proteinDescription', 'genes', 'comments', 'features', 'keywords',
    # 'references', 'uniProtKBCrossReferences', 'sequence', 'extraAttributes'])

    # output_info = np.load('dict.npy', allow_pickle=True).item()
    # pprint.pprint(output_info)
