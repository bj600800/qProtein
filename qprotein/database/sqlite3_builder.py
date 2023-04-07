"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: Base class for SQL builder class
# ------------------------------------------------------------------------------
"""
import sys
sys.path.insert(0, '../')
from abc import abstractmethod
import sqlite3
import gzip
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


class SqlBuilder(object):

    def __init__(self, sql_db, table_name, column_definition):
        self.sql_path = sql_db
        self.table_name = table_name
        self.column_definition = column_definition

        self.records = []
        self.record_counter = 0
        self.total_records = 0

    @staticmethod
    def read_text_generator(filename, header):
        with open(filename, 'r') as f:
            if header == True or 'T':
                next(f)
            for line in f:
                yield line

    @staticmethod
    def read_gz_generator(file):
        with gzip.open(file, 'rt') as uniprot_dat:
            for line in uniprot_dat:
                yield line

    @staticmethod
    def create_table(table_name, sql_db, sql):
        connect = sqlite3.connect(sql_db)
        cursor = connect.cursor()
        cursor.execute(f"DROP TABLE IF EXISTS {table_name}")
        cursor.execute(sql)
        return cursor

    def add_column(self, cursor, table_name, column_definition: list):
        cursor.execute(f"SELECT * FROM {table_name}")
        existing_columns = [row[0] for row in cursor.description]
        new_columns = [col for col in column_definition if col[0] not in existing_columns]
        if not new_columns:
            return
        for name, dtype in new_columns:
            sql = f"ALTER TABLE {table_name} ADD COLUMN {name} {dtype}"
            cursor.execute(sql)

    @staticmethod
    def connect_sql(sql_db):
        connect = sqlite3.connect(sql_db)
        cursor = connect.cursor()
        return cursor

    def insert_many(self, cursor, records: list, values_num):
        cursor.execute("begin")
        question_mark = ', '.join(['?'] * values_num)
        cursor.executemany(f"INSERT INTO {self.table_name} VALUES({question_mark})", records)
        cursor.execute("commit")

    def insert_many_columns(self, cursor, columns, records: list, values_num):
        columns = ', '.join(columns)
        cursor.execute("begin")
        question_mark = ', '.join(['?'] * values_num)
        cursor.executemany(f"INSERT INTO {self.table_name} ({columns}) VALUES ({question_mark})", records)
        cursor.execute("commit")

    def update_many_columns(self, cursor, table_name, columns, records: list):
        cursor.execute("begin")
        target_column = columns[0]
        columns_question = ', '.join([i+'=?' for i in columns[1:]])
        format_records = [i[1:] + (i[0],) for i in records]
        sql = f"UPDATE {table_name} SET {columns_question} where {target_column}=?"
        cursor.executemany(sql, format_records)
        cursor.execute("commit")

    def insert_update_columns(self, cursor, table_name, columns, records: list):
        batch_size = 500000
        total_query = len(records)
        question_num = len(columns)
        column_name = ', '.join(columns)
        for i in range(0, total_query, batch_size):
            batch_acc = records[i:i+batch_size]
            question_mark = ', '.join(['?'] * question_num)
            excluded_clause = ', '.join([i+'=excluded.'+i for i in columns[1:]])

            # use this sql need has the primary key.
            sql_cmd = f"INSERT INTO {table_name} ({column_name}) VALUES ({question_mark})" \
                      f"ON CONFLICT ({columns[0]}) DO UPDATE SET {excluded_clause}"
            cursor.execute('begin')
            cursor.executemany(sql_cmd, batch_acc)
            cursor.execute('commit')

    def create_index(self, cursor, columns):
        for column in columns:
            cursor.execute('begin')
            cursor.execute(f"CREATE INDEX IF NOT EXISTS index_{column} ON {self.table_name} ({column})")
            cursor.execute('commit')
        cursor.close()

    @abstractmethod
    def run(self):
        """
        This is an abstract method that needs to be implemented in the subclass.
        """
        raise NotImplementedError("The read method must be implemented by a subclass")
