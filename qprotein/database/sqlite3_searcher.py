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