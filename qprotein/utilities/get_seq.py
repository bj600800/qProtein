"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/14

# Description: 
# ------------------------------------------------------------------------------
"""
import os

from qprotein.database.sqlite3_searcher import SqlSearch


def get_query(cif_dir):
    query_list = [os.path.splitext(i)[0].split("#")[0] for i in os.listdir(cif_dir)]
    return query_list


def get_seq(query_list, fasta_sql, task_name):
    cursor = SqlSearch.connect_sql(fasta_sql)
    sql_cmd = "SELECT query_name, sequence FROM {} WHERE query_name IN ({})" \
        .format(task_name, ', '.join(['?'] * len(query_list)))
    fasta = SqlSearch.fetch_many_results(cursor, sql_cmd, query_list)
    return fasta


def write_fasta(fasta, fasta_file):
    with open(fasta_file, 'w') as wf:
        for name, seq in tqdm(fasta):
            wf.write('>' + name + '\n')
            wf.write(seq + '\n')


if __name__ == '__main__':
    root_dir = r"D:\subject\active\1-qProtein\data"
    task_name = ["tibet", "manure"]
    n = 0
    work_dir = os.path.join(root_dir, task_name[n])
    structure_dir_path = os.path.join(work_dir, 'high_plddt')
    fasta_sql = os.path.join(root_dir, 'qprotein_db.db')
    fasta_file = os.path.join(work_dir, task_name[n]+'_high_plddt.fasta')
    query_list = get_query(structure_dir_path)
    fasta = get_seq(query_list, fasta_sql, task_name[n])
    write_fasta(fasta, fasta_file)
