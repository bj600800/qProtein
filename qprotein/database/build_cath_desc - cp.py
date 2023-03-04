# -*- coding: utf-8 -*-
# @Time    : 2023/2/24 11:34
# @Author  : Zhixin Dou
"""
Build cath description sqlite file
"""
import re
import sqlite3
import time

from qprotein.utilities import log

logger = log.setup_log(name='sql_Builder')


class SqlBuilder(object):

    def __init__(self, desc_file, sql_db, table_name, column_definition):
        self.desc_file = desc_file
        self.sql_db = sql_db
        self.table_name = table_name
        self.column_definition = column_definition
        self.records = []
        self.main()

    def read_desc_file(self):
        logger.info(f"Step 1 -> Read source file: {self.desc_file}")

        with open(self.desc_file, 'r') as f:
            desc_content = [line for line in f.readlines() if not line.startswith('#')]
        return desc_content

    def create_cath_sql(self):
        logger.info(f"Step 2 -> Create SQL table: {self.table_name}")

        connect = sqlite3.connect(self.sql_db)
        cursor = connect.cursor()
        cursor.execute(f"DROP TABLE IF EXISTS {self.table_name}")
        definition_str = ','.join([' '.join(i) for i in self.column_definition.items()])
        cursor.execute(f"CREATE TABLE {self.table_name} ({definition_str})")
        return cursor

    def pack_records(self, desc_content):
        # TODO: overwrite in subclass
        # Different object has different data structure.
        logger.info("Step 3 -> Pack records for SQL batch inserting")
        record_items = {}
        prefixes = ["DOMAIN", "CATHCODE", "CLASS", "ARCH", "TOPOL", "HOMOL"]
        reg = re.compile(r"  +")

        for line in desc_content:
            for prefix in prefixes:
                if line.startswith(prefix):
                    record_items[prefix] = reg.split(line)[1].rstrip()
            if line.startswith("//\n"):
                self.records.append(tuple(record_items.values()))

    def insert_batch(self, cursor, records: list):
        logger.info(f"Step 4 -> Insert data in bulk to SQL table: {self.table_name}")
        cursor.execute("begin")
        cursor.executemany(f"INSERT INTO {self.table_name} VALUES(?, ?, ?, ?, ?, ?)", records)
        cursor.execute("commit")

    def create_index(self, cursor):
        logger.info(f'Step 5 -> Creating index for SQL table: {self.table_name}')
        column_name = ', '.join(self.column_definition.keys())
        cursor.execute(f"CREATE INDEX cath_desc_idx ON cath_desc ({column_name})")

    def exec_sql(self, desc_content, cursor):
        self.pack_records(desc_content)
        self.insert_batch(cursor, self.records)
        self.create_index(cursor)

    def main(self):
        logger.info(f'Start to build SQL database for {table_name}')
        desc_content = self.read_desc_file()
        cursor = self.create_cath_sql()
        self.exec_sql(desc_content, cursor)
        logger.info(f'Successfully built SQL database for {table_name}')


if __name__ == '__main__':
    desc_file = r'G:\DB\cath\cath-classification-data_txt\cath-domain-description-file.txt'
    sql_db = r'C:\Users\bj600\Desktop\cath.db'
    table_name = 'cath_desc'
    column_definition = {"domain_name": "TEXT", "cath_code": "TEXT",
                         "cath_class": "TEXT", "cath_architecture": "TEXT",
                         "cath_topology": "TEXT", "cath_homologous": "TEXT"
                         }
    cath_builder = SqlBuilder(desc_file, sql_db, table_name, column_definition)
