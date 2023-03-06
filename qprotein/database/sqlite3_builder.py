"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: Base class for SQL builder class
# ------------------------------------------------------------------------------
"""

from abc import abstractmethod
import sqlite3
import gzip


class SqlBuilder(object):

    def __init__(self, sql_db, table_name, column_definition):
        self.sql_db = sql_db
        self.table_name = table_name
        self.column_definition = column_definition

        self.records = []
        self.record_counter = 0

    @staticmethod
    def read_text_generator(filename):
        with open(filename, 'r') as f:
            for line in f:
                yield line

    @staticmethod
    def read_gz_generator(file):
        with gzip.open(file, 'rt') as uniprot_dat:
            for line in uniprot_dat:
                yield line

    def create_table(self):
        connect = sqlite3.connect(self.sql_db)
        cursor = connect.cursor()
        cursor.execute(f"DROP TABLE IF EXISTS {self.table_name}")
        definition_str = ','.join([' '.join(i) for i in self.column_definition.items()])
        cursor.execute(f"CREATE TABLE {self.table_name} ({definition_str})")
        return cursor

    def insert_many(self, cursor, records: list, values_num):
        cursor.execute("begin")
        question_mark = ', '.join(['?'] * values_num)
        cursor.executemany(f"INSERT INTO {self.table_name} VALUES({question_mark})", records)
        cursor.execute("commit")

    def create_index(self, cursor):
        column_name = ', '.join(self.column_definition.keys())
        cursor.execute(f"CREATE INDEX {self.table_name}_idx ON {self.table_name} ({column_name})")

    @abstractmethod
    def run(self):
        """
        This is an abstract method that needs to be implemented in the subclass.
        """
        raise NotImplementedError("The read method must be implemented by a subclass")
