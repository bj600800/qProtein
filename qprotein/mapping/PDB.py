"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/07/07

# Description: Step 1. Get structures from Protein Data Bank.
# ------------------------------------------------------------------------------
"""
import os
import warnings
import tempfile
import numpy as np
import requests
from requests.adapters import Retry
import biotite.structure as struc
import biotite.structure.io as strucio
from Bio.SeqUtils import IUPACData
from qprotein.preprocessing.sqlite3_searcher import SqlSearch
from qprotein.utilities import logger
from Bio import BiopythonWarning
warnings.simplefilter("ignore", BiopythonWarning)
logger = logger.setup_log(name=__name__)


def start_request_session():
    retries = Retry(total=3, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504])
    session = requests.Session()
    session.mount('http://', requests.adapters.HTTPAdapter(max_retries=retries))
    session.mount('https://', requests.adapters.HTTPAdapter(max_retries=retries))
    return session


class PDBMapper(SqlSearch):
    def __init__(self, sql_db, folder):
        self.sql_db = sql_db
        self.positive_folder = os.path.join(folder, "positive")
        self.negative_folder = os.path.join(folder, "negative")
        if not os.path.exists(self.positive_folder):
            os.mkdir(self.positive_folder)
        if not os.path.exists(self.negative_folder):
            os.mkdir(self.negative_folder)
            
    def get_pdb_id(self):
        cursor = self.connect_sql(self.sql_db)
        sql_cmd = "SELECT query_name, accession_id, subject_start, subject_end, label From results_summary WHERE source == 'PDB'"
        pdb_ret = self.fetch_results(cursor, sql_cmd)
        return pdb_ret

    def get_exist_query(self):
        exist_structure = [file for folder in [self.positive_folder, self.negative_folder] for file in os.listdir(folder) if
                           os.path.isfile(os.path.join(folder, file)) and os.path.getsize(
                                   os.path.join(folder, file)
                                   ) > 0]
        exist_query_name = [os.path.splitext(item)[0] for item in exist_structure]
        return exist_query_name
    
    @staticmethod
    def get_missing_residue_number(text, chain_id, subject_start):
        aa_list = [value.upper() for value in IUPACData.protein_letters_1to3.values()]
    
        # 查找并处理MISSING RESIDUES信息
        missing_residues = []
        for line in text.split('\n'):
            if line.startswith('REMARK 465') and len(line.split()) == 5:
                residue = line.split()[2]
                chain = line.split()[3]
                if residue in aa_list and chain == chain_id:
                    missing = line.split()[2:][::2]
                    missing[0] = IUPACData.protein_letters_3to1[missing[0][0] + missing[0][1:].lower()]
                    if int(missing[1]) < subject_start:
                        missing_residues.append(missing)
        return len(missing_residues)
        
    def get_structure(self, query_info):
        query_name = query_info[0]
        pdb_full_id = query_info[1]
        start_residue = query_info[2]
        end_residue = query_info[3]
        label = query_info[4]
        pdb_id, chain_id = pdb_full_id.split("_")
        # 构建PDB的API链接
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/70.0.3538.25 Safari/537.36 Core/1.70.3861.400 QQBrowser/10.7.4313.400",
            "From": "bj600800@gmail.com"  # ALLWAYS TELL WHO YOU ARE
            }
        api_url = f"https://files.rcsb.org/view/{pdb_id}.pdb"
        session = start_request_session()
        # 发送GET请求并获取结构数据
        try:
            response = session.get(api_url, headers=headers)
            pdb_string = response.text
    
            # 将结构数据保存为临时PDB文件
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".pdb") as temp:
                temp.write(pdb_string)
                temp_pdb_filepath = temp.name
                temp.flush()
    
            # 解析PDB结构数据
            structure = strucio.load_structure(temp_pdb_filepath)
            
            # remove non-amino acids
            mask_residue = struc.filter_amino_acids(structure)
            structure = structure[mask_residue]
            
            # get the target chain
            mask_chain = (structure.chain_id == chain_id)
            structure = structure[mask_chain]
            subject_start = min(structure.res_id)
            
            # get missing num
            missing_num = self.get_missing_residue_number(pdb_string, chain_id=chain_id, subject_start=subject_start)
            start_residue -= missing_num
            end_residue -= missing_num
            
            # renumber residue id for wrong ids
            structure = struc.renumber_res_ids(structure, start=1)
            
            # make path
            if label == "positive":
                clean_file_path = os.path.join(self.positive_folder, f"{query_name}.pdb")
            else:
                clean_file_path = os.path.join(self.negative_folder, f"{query_name}.pdb")

            # filter given start and end position residues
            mask_segment = np.isin(structure.res_id, list(range(start_residue, end_residue+1)))
            structure = structure[mask_segment]
            
            # renumber residue ids after segmenting structure
            structure = struc.renumber_res_ids(structure, start=1)
            
            # save structure
            strucio.save_structure(clean_file_path, structure)
            os.remove(temp_pdb_filepath)
            pdb_info = [query_name, clean_file_path]
            print("Structure path: ", clean_file_path)
            # return pdb_info
        except:
            pass
        
    def run(self):
        query_info = []
        pdb_ret = self.get_pdb_id()
        exist_query_name = self.get_exist_query()
        pdb_query = [result for result in pdb_ret if result[0] not in exist_query_name]
        for one_ret in pdb_query:
            query_info.append(self.get_structure(one_ret))
        # return query_info
    
if __name__ == '__main__':
    structure_folder = r"D:\subject\active\1-qProtein\data\enzymes\GH10\2_StrucMapping"
    summary_sql_path = r"D:\subject\active\1-qProtein\data\enzymes\GH10\1_preprocessing\qprotein_results.db"
    if not os.path.exists(structure_folder):
        os.mkdir(structure_folder)
    test = PDBMapper(sql_db=summary_sql_path, folder=structure_folder)
    test.run()
    