"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: Base class for SQL searcher class
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

    @staticmethod
    def connect_sql(sql_db):
        connect = sqlite3.connect(sql_db)
        cursor = connect.cursor()
        return cursor

    @staticmethod
    def fetch_results(cursor, sql_cmd):
        cursor.execute(sql_cmd)
        results_list = cursor.fetchall()
        return results_list

    @staticmethod
    def fetch_many_results(cursor, sql_cmd, search_targets):
        cursor.execute(sql_cmd, search_targets)
        results_list = cursor.fetchall()
        return results_list

    # def update_many(self, cursor, table_name, set_column_name, condition_column_name, data):
    #     sql_cmd = f'UPDATE {table_name} SET {set_column_name}=? WHERE {condition_column_name}=?'
    #     cursor.executemany(sql_cmd, data)
    #     cursor.execute("commit")

    @staticmethod
    def read_text_generator(filename):
        with open(filename, 'r') as f:
            for line in f:
                yield line


class AFmetaSearch(SqlSearch):
    def __init__(self, sql_db, table_name):
        super(AFmetaSearch, self).__init__(sql_db)

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

# if __name__ == '__main__':
#     table_name = 'AFmeta'
#     AF_sql = r'G:\DB\AFmeta.db'
#     AFsearch = AFmetaSearch(AF_sql, table_name)
#     AFsearch.run()
