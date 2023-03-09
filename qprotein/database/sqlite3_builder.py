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

from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)

class SqlBuilder(object):

    def __init__(self, sql_db, table_name, column_definition):
        self.sql_db = sql_db
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
    def create_table(table_name, sql_db, column_definition):
        connect = sqlite3.connect(sql_db)
        cursor = connect.cursor()
        cursor.execute(f"DROP TABLE IF EXISTS {table_name}")
        definition_str = ', '.join([' '.join(i) for i in column_definition])
        cursor.execute(f"CREATE TABLE {table_name} ({definition_str})")
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

    def update_many_columns(self, cursor, table_name, columns, records: list, target_column, values_num):
        # TODO the method should be separated into Update and insert

        update_records = []
        insert_records = []
        search_ret = []
        batch_size = 100000
        target_value = [record[0] for record in records]
        total_query = len(target_value)
        for i in range(0, total_query, batch_size):
            batch_acc = target_value[i:i+batch_size]
            question_mark = ', '.join(['?'] * len(batch_acc))
            sql_cmd = f"SELECT {target_column} FROM {table_name} WHERE {target_column} IN ({question_mark})"
            cursor.execute(sql_cmd, batch_acc)
            ret = cursor.fetchall()
            search_ret.extend([str(i[0]) for i in ret])

        for record in records:
            if record[0] in search_ret:
                update_records.append(record)
            else:
                insert_records.append(record)

        cursor.execute("begin")
        if update_records:
            logger.info('Has {} records to update.'.format(len(update_records)))
            update_records = [(t[1:] + (t[0],)) for t in update_records]
            columns_question = ', '.join([i+' = ?' for i in columns[1:]])
            cursor.executemany(f"UPDATE {self.table_name} SET {columns_question} WHERE {columns[0]} = ?", update_records)

        if insert_records:
            logger.info('Has {} records to insert.'.format(len(insert_records)))
            columns = ', '.join(columns)
            question_mark = ', '.join(['?'] * values_num)
            cursor.executemany(f"INSERT INTO {self.table_name} ({columns}) VALUES ({question_mark})", insert_records)
        cursor.execute("commit")

    def create_index(self, cursor, idx_name, column_name):
        sql = f"CREATE INDEX {idx_name}_idx ON {self.table_name} ({column_name})"
        cursor.execute(sql)


    @abstractmethod
    def run(self):
        """
        This is an abstract method that needs to be implemented in the subclass.
        """
        raise NotImplementedError("The read method must be implemented by a subclass")
