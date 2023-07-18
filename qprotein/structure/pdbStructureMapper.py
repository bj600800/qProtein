"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/07/07

# Description: 
# ------------------------------------------------------------------------------
"""
import os
import warnings
import requests
from Bio.PDB import PDBParser, PDBIO, Select, standard_aa_names


from qprotein.database.sqlite3_builder import SqlBuilder
from qprotein.database.sqlite3_searcher import SqlSearch
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)
from Bio import BiopythonWarning
warnings.simplefilter("ignore", BiopythonWarning)


# 清洗选择器，去除水分子、小分子和杂原子，只保留指定链的氨基酸残基
class CleanSelector(Select):
    def __init__(self, chain_id, start_residue, end_residue):
        self.chain_id = chain_id
        self.start_residue = start_residue
        self.end_residue = end_residue
    
    def accept_model(self, model):
        return True
    
    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id
    
    def accept_residue(self, residue):
        if residue.get_parent().get_id() != self.chain_id:
            return False
        if residue.resname in standard_aa_names:
            residue_id = int(residue.get_id()[1])
            if self.start_residue <= residue_id <= self.end_residue:
                return True
        return False

def get_pdb_id(blast_out):
    all_pdb_parm = []
    with open(blast_out, "r") as f:
        content = f.readlines()
    for line in content:
        all_pdb_parm.append([line.split("\t")[0], line.split("\t")[1], [int(line.split("\t")[7]), int(line.split("\t")[8])]])
    return all_pdb_parm

def get_complete_sequence_from_pdb(text):
    """
    ERROR: 对于PDB缺失部分残基导致无法根据序列比对结果准确截取结构域的情况，暂时搁置。
    """
    # 查找并处理MISSING RESIDUES信息
    missing_residues = []
    for line in text.split('\n'):
        if line.startswith('REMARK 465') and 'MISSING RESIDUES' in line:
            residue_info = line.split()[3:]
            for i in range(0, len(residue_info), 2):
                residue_id = residue_info[i + 1]
                chain_id = residue_info[i][0]
                missing_residues.append((residue_id, chain_id))
    
    return missing_residues


def adjust_residue_number(structure):
    """
    调整开始序数大于1的情况
    """
    residues = list(structure.get_residues())
    start_resseq = [residue.get_id()[1] for residue in residues][0]
    if start_resseq > 1:
        diff = 1 - start_resseq
        for residue in residues:
            parent_chain = residue.get_parent()
            parent_chain.detach_child(residue.get_id())
            new_id = (residue.id[0], residue.id[1]+diff, residue.id[2])
            residue.id = new_id
            parent_chain.add(residue)
    return structure
    
    
def get_structure(pdb_parm, folder):
    query_name = pdb_parm[0]
    pdb_full_id = pdb_parm[1]
    start_residue, end_residue = pdb_parm[2]
    pdb_id, chain_id = pdb_full_id.split("_")
    
    # 构建PDB的API链接
    api_url = f"https://files.rcsb.org/view/{pdb_id}.pdb"
    
    # 发送GET请求并获取结构数据
    response = requests.get(api_url)
    structure_data = response.text
    
    # 将结构数据保存为临时PDB文件
    temp_pdb_filename = f"{pdb_id}.pdb"
    with open(temp_pdb_filename, "w") as f:
        f.write(structure_data)
    
    # 创建PDB解析器
    parser = PDBParser()
    
    # 解析PDB结构数据
    structure = parser.get_structure(pdb_id, temp_pdb_filename)

    # 调整残基的序号
    structure = adjust_residue_number(structure)

    # 使用选择器筛选结构中的链和原子，并保存为新的PDB文件
    if not os.path.exists(folder):
        os.makedirs(folder)
    clean_pdb_file = os.path.join(folder, f"{query_name}&PDB&{pdb_full_id}.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save(clean_pdb_file, CleanSelector(chain_id=chain_id, start_residue=start_residue, end_residue=end_residue))
    os.remove(temp_pdb_filename)
    pdb_info = [query_name, clean_pdb_file]
    # 打印结果
    print(f"PDB structure with {pdb_full_id} has been downloaded and cleaned.")
    print(f"The cleaned structure is saved as {clean_pdb_file}.")
    return pdb_info


class PDBUpdate(SqlBuilder, SqlSearch):
    def __init__(self, summary_sql_path, pdb_info_list):
        super().__init__()
        self.summary_sql_path = summary_sql_path
        self.pdb_info_list = pdb_info_list
        self.table_name = "results_summary"
        self.column_definition = [('query_name', 'TEXT'), ('pdb_file', 'TEXT')]
        self.record_items = {column_name[0]: '' for column_name in self.column_definition}
        self.records = []
        self.record_counter = 0
        self.total_records = 0

    def parse_pdbInfo(self, pdb_info):
        self.record_items['query_name'] = pdb_info[0]
        self.record_items['pdb_file'] = pdb_info[1]
        
    def format_pdbInfo(self, pdb_info):
        # prepare
        self.parse_pdbInfo(pdb_info)
        self.records.append(tuple(self.record_items.values()))
        self.record_counter += 1
        self.total_records += 1

    def operate_pdbInfo(self, cursor, columns):
        for pdb_info in self.pdb_info_list:
            if self.record_counter == 500000:
                self.update_many_columns(
                    cursor=cursor, table_name="results_summary",
                    columns=columns, records=self.records
                    )
                logger.info(f"Insert 500000 records. Total records: {str(self.total_records)}")
                self.records = []
                self.record_counter = 0
            else:
                self.format_pdbInfo(pdb_info)
    
        if self.record_counter > 0:
            logger.info(f"Insert {self.record_counter} records")
            print(self.records)
            self.update_many_columns(cursor=cursor, table_name="results_summary",
                                     columns=columns, records=self.records)
        logger.info(f"Total records:" + str(self.total_records))

    def run(self):
        logger.info(f"Update pdb structure infomation to summary database")
        columns = [i[0] for i in self.column_definition]
        cursor = self.connect_sql(sql_db=self.summary_sql_path)
        self.add_column(cursor=cursor, table_name="results_summary", column_definition=self.column_definition)
        self.operate_pdbInfo(cursor=cursor, columns=columns)
        self.create_index(cursor=cursor, columns=columns, table_name=self.table_name)
        
def main():
    work_dir = r"D:\subject\active\1-qProtein\data\enzymes\cazy_endo-1_4-beta-xylanase"
    summary_sql_path = r"D:\subject\active\1-qProtein\data\enzymes\cazy_endo-1_4-beta-xylanase\qprotein_results.db"
    blast_out = os.path.join(work_dir, "11_map_pdb_sequence_filtered.out")
    structure_folder = os.path.join(work_dir, "structure")
    all_pdb_parm = get_pdb_id(blast_out=blast_out)
    pdb_info_list = []
    for pdb_parm in all_pdb_parm:
        pdb_info = get_structure(pdb_parm=pdb_parm, folder=structure_folder)
        pdb_info_list.append(pdb_info)
    pdbUpdater = PDBUpdate(summary_sql_path=summary_sql_path, pdb_info_list=pdb_info_list)
    pdbUpdater.run()
    
if __name__ == '__main__':
    main()
    