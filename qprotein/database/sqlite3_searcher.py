"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: Base class for SQL builder class, subclass AFmeta, uniprot
# ------------------------------------------------------------------------------
"""

from abc import abstractmethod
import sqlite3

from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


class SqlSearch(object):
    def __init__(self, sql_db):
        self.records = []
        self.record_counter = 0
        self.sql_db = sql_db

    def connect_sql(self, sql_db):
        connect = sqlite3.connect(sql_db)
        cursor = connect.cursor()
        return cursor

    def fetch_results(self, cursor, sql_cmd):
        ret = cursor.execute(sql_cmd)
        results_list = ret.fetchall()
        return results_list

    @staticmethod
    def read_text_generator(filename):
        with open(filename, 'r') as f:
            for line in f:
                yield line


class AFmetaSearch(SqlSearch):
    def __init__(self, sql_db, table_name):
        super().__init__(sql_db)

        self.sql_db = sql_db
        self.table_name = table_name

    def search_sql(self, cursor, target_list):
        results = []
        for i in target_list:
            sql_query = f"SELECT * FROM {self.table_name} WHERE UniDesc=? Limit 1"
            cursor.execute(sql_query, (i[0],))
            _results = cursor.fetchall()
            results.append(_results)
        return results

    def write_output(self, results, file):
        with open(file, 'w') as wt:
            for result in results:
                for one_result in result:
                    str_result = [str(i) for i in one_result]
                    wt.write('\t'.join(str_result) + '\n')

    def run(self):
        cursor = self.connect_sql(self.sql_db)
        sql_target = f"SELECT UniDesc, COUNT (UniDesc) FROM {self.table_name} GROUP BY UniDesc ORDER BY COUNT (UniDesc) DESC LIMIT 26"
        logger.info(f"Fetch targets from {self.table_name}")
        target_list = self.fetch_results(cursor, sql_target)

        logger.info("Search results")
        results = self.search_sql(cursor, target_list)
        file = r'D:\subject\active\1-qProtein\data\alphafold_metadata\two-results_for_function_cluster.txt'

        logger.info(f"Write the results in {file}")
        self.write_output(results, file)
        cursor.close()


# if __name__ == '__main__':
#     table_name = 'AFmeta'
#     AF_sql = r'G:\DB\AFmeta.db'
#     AFsearch = AFmetaSearch(AF_sql, table_name)
#     AFsearch.run()


class UniprotSearch(SqlSearch):
    def __init__(self, sql_db, results_db):
        super().__init__(sql_db)

        self.results_db = results_db

    def fetch_accession(self, results_db, table_name):
        cursor_results = self.connect_sql(results_db)
        sql_cmd = f"SELECT sprot100_acc FROM {table_name} WHERE sprot100_acc != ''"
        accession_list = [i[0] for i in self.fetch_results(cursor_results, sql_cmd)]
        return accession_list

    def get_annotation(self, sql_db, table_name, accession_list):
        cursor_db = self.connect_sql(sql_db)
        acc_id = ','.join(['?'] * len(accession_list))
        cursor_db.execute(f"SELECT * FROM {table_name} WHERE accession in ({acc_id})", accession_list)
        annotation = cursor_db.fetchall()
        return annotation

    def run(self):
        accession_list = self.fetch_accession(self.results_db, table_name='summary')
        annotation = self.get_annotation(self.sql_db, 'uniprot_sprot', accession_list)
        print(annotation)
        return annotation

if __name__ == '__main__':
    results_db = r'D:\subject\active\1-qProtein\data\tibet\summary.db'
    sql_db = r'G:\DB\uniprot_sprot.db'
    sprot_dat = UniprotSearch(sql_db, results_db)
    sprot_dat.run()
