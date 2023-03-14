"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/12

# Description: Fetch target sequences from query_seq database with or without qprotein annotation.
# ------------------------------------------------------------------------------
"""

from sqlite3_searcher import SqlSearch


def fetch_target(sql_db, column, key_word):
    """
    return: {unavailable: {query_name__family: sequence}
             available: {query_name__family: {sprot_acc: sequence}, {trembl_acc: sequence}}}
    """
    cursor = SqlSearch.connect_sql(sql_db)
    sql = f"SELECT query_name FROM results_summary WHERE {column} LIKE '%{key_word}%'"
    target_name = [i[0] for i in SqlSearch.fetch_results(cursor, sql)]

    return target_name


def fetch_sequence_of_target(sql_fasta, target_name):
    """
    Prepare sequences for domain annotation and segmentation with qprotein annotation
    """
    cursor = SqlSearch.connect_sql(sql_fasta)
    batch_size = 100000
    target_protein = []
    total_query = len(target_name)
    for i in range(0, total_query, batch_size):
        batch_acc = target_name[i:i + batch_size]
        sql_cmd = "SELECT query_name, sequence, length(sequence) FROM query_seq " \
                  "WHERE query_name IN ({}) AND length(sequence) > 200"\
                  .format(', '.join(['?'] * len(batch_acc)))
        cursor.execute(sql_cmd, batch_acc)
        return_dat = cursor.fetchall()
        target_protein.extend(return_dat)

    return target_protein


def fetch_sequence(sql_fasta):
    """
    without qprotein annotation
    """
    pass



if __name__ == '__main__':
    sql_db = r'D:\subject\active\1-qProtein\data\tibet\qprotein_results.db'
    sql_fasta = r'D:\subject\active\1-qProtein\data\tibet\qprotein_db.db'
    column_key = [('cazy_family', 'GH'), ('cazy_family', 'GT'), ('cazy_family', 'AA'),
                  ('cazy_family', 'PL'), ('cazy_family', 'CE'), ('merops_family', '')]
    protein = []
    for column, key_word in column_key:
        result_list = fetch_target(sql_db, column, key_word)
        target_protein = fetch_sequence_of_target(sql_fasta, result_list)
        protein.extend(target_protein)

    with open(r'D:\subject\active\1-qProtein\data\tibet\tibet_target_ident90.fasta', 'w') as wf:
        for i in protein:
            wf.write('>'+i[0]+'\n')
            wf.write(i[1]+'\n')
