# -*- coding: utf-8 -*-
# @Time    : 2023/2/24 11:34
# @Author  : Zhixin Dou
"""
Build cath description sqlite file
"""
import re

from qprotein.utilities import logger
from qprotein.database.sqlite3_builder import SqlBuilder

logger = logger.setup_log(name=__name__)
print("-" * 60)
print("The following logs from the logger: ", logger.name)


class CATHsql(SqlBuilder):

    def __init__(self, file, sql_db, table_name, column_definition):
        super().__init__(sql_db, table_name, column_definition)

        self.file = file
        self.records = []

    def read_file(self):
        with open(self.file, 'r') as f:
            desc_content = [line for line in f.readlines() if not line.startswith('#')]
        return desc_content

    def pack_records(self, desc_content):
        # TODO: overwrite in subclass
        # Different object has different data structure.
        record_items = {}
        prefixes = ["DOMAIN", "CATHCODE", "CLASS", "ARCH", "TOPOL", "HOMOL"]
        reg = re.compile(r"  +")

        for line in desc_content:
            for prefix in prefixes:
                if line.startswith(prefix):
                    record_items[prefix] = reg.split(line)[1].rstrip()
            if line.startswith("//\n"):
                self.records.append(tuple(record_items.values()))

    def run(self):
        logger.info(f"Start to create SQL table: {table_name} in SQL file {sql_db}")

        logger.info(f"Step 1 -> Read source file: {self.file}")
        desc_content = self.read_file()

        logger.info(f"Step 2 -> Create SQL table: {table_name}")
        cursor = self.create_cath_sql()

        logger.info("Step 3 -> Pack records for SQL batch inserting")
        self.pack_records(desc_content)

        logger.info(f"Step 4 -> Insert data in bulk to SQL table: {table_name}")
        self.insert_batch(cursor, self.records)

        logger.info(f'Step 5 -> Creating index for SQL table: {table_name}')
        self.create_index(cursor)

        logger.info(f'Successfully built SQL database for {table_name}')


if __name__ == '__main__':
    desc_file = r'G:\DB\cath\cath-classification-data_txt\cath-domain-description-file.txt'
    sql_db = r'C:\Users\bj600\Desktop\cath.db'
    table_name = 'cath_desc'
    column_definition = {"domain_name": "TEXT", "cath_code": "TEXT",
                         "cath_class": "TEXT", "cath_architecture": "TEXT",
                         "cath_topology": "TEXT", "cath_homologous": "TEXT"
                         }
    cath_builder = CATHsql(desc_file, sql_db, table_name, column_definition)
    cath_builder.run()
