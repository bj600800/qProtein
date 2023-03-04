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

        self.build_sqlite()

    def create_cath_sql(self):
        connect = sqlite3.connect(self.sql_db)
        cursor = connect.cursor()
        cursor.execute(f"DROP TABLE IF EXISTS {self.table_name}")
        cursor.execute(f"CREATE TABLE {self.table_name} ({self.column_definition})")
        return cursor

    def read_desc_file(self):
        with open(self.desc_file, 'r') as f:
            desc = [line for line in f.readlines() if not line.startswith('#')]
        return desc

    def build_sqlite(self):
        cursor = self.create_cath_sql()
        records = []
        record_counter = 0
        column_items = {}
        reg = re.compile(r"  +")
        prefixes = ["DOMAIN", "CATHCODE", "CLASS", "ARCH", "TOPOL", "HOMOL"]

        for line in self.read_desc_file():
            if record_counter == 50000:
                print('Insert 50000 items.')
                cursor.execute("begin")
                cursor.executemany(f"INSERT INTO {self.table_name} VALUES(?, ?, ?, ?, ?, ?)", records)
                cursor.execute("commit")
                records = []
                record_counter = 0
                records.append(tuple(column_items.values()))
                record_counter += 1
            else:
                for prefix in prefixes:
                    if line.startswith(prefix):
                        column_items[prefix] = reg.split(line)[1].rstrip()
                if line.startswith("//\n"):
                    records.append(tuple(column_items.values()))
                    record_counter += 1

        if record_counter > 0:
            print(records)
            print('Procession remaining records.')
            cursor.execute("begin")
            cursor.executemany(f"INSERT INTO {self.table_name} VALUES(?, ?, ?, ?, ?, ?)", records)
            cursor.execute("commit")
        print(f'Creating index for {self.table_name}.')
        cursor.execute(f"CREATE INDEX cath_desc_idx ON cath_desc (domain_name, cath_code, cath_class,"
                       "cath_architecture, cath_topology, cath_homologous)"
                       )

        print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
        print('Finished sql_db building for cath_desc.')


if __name__ == '__main__':
    desc_file = r'G:\DB\cath\cath-classification-data_txt\cath-domain-description-file.txt'
    sql_db = r'G:\DB\cath\cath-classification-data_txt\cath.db'
    table_name = 'cath_desc'
    column_definition = 'domain_name TEXT, cath_code TEXT, cath_class TEXT,' \
                        'cath_architecture TEXT, cath_topology TEXT, cath_homologous TEXT'
    cath_builder = SqlBuilder(desc_file, sql_db, table_name, column_definition)
