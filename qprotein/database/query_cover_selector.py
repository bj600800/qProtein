"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/04/07

# Description: 
# ------------------------------------------------------------------------------
"""
from qprotein.database.sqlite3_searcher import SqlSearch


def get_name_match_length(trembl_out):
    with open(trembl_out, 'r') as f:
        content = f.readlines()
    query_match = [(i.split('\t')[0], i.split('\t')[3]) for i in content]  # name, match_length
    return query_match


def get_fasta(fasta_sql, query_match):
    cursor = SqlSearch.connect_sql(fasta_sql)
    query_name = [i[0] for i in query_match]
    sql_cmd = "SELECT query_name, sequence from query_seq where query_name in ({})" .format(', '.join(['?'] * len(query_name)))
    results = SqlSearch.fetch_many_results(cursor, sql_cmd, query_name)
    fasta = [(result[0], result[1]) for result in results for query in query_match
             if result[0] == query[0] and len(result[1]) > 200 and int(query[1]) / len(result[1]) >= 0.9]
    return fasta


def write_fasta(fasta, write_path):
    with open(write_path, 'w') as wf:
        for i in fasta:
            wf.write('>' + i[0] + '\n')
            wf.write(i[1] + '\n')


def run():
    trembl_out = r'D:\subject\active\1-qProtein\data\tibet\trembl_90.tab'
    fasta_sql = r'D:\subject\active\1-qProtein\data\qprotein_db.db'
    write_path = r'D:\subject\active\1-qProtein\data\tibet\trembl90_90.fasta'
    query_match = get_name_match_length(trembl_out)
    fasta = get_fasta(fasta_sql, query_match)
    write_fasta(fasta, write_path)


run()
