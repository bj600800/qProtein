# -*- coding: utf-8 -*-
# @Time    : 2023/2/24 11:34
# @Author  : Zhixin Dou
"""
Build cath description sqlite file, source file read in by a generator from SqlBuilder
"""
import re

from qprotein.utilities import logger
from qprotein.database.sqlite3_builder import SqlBuilder

logger = logger.setup_log(name=__name__)


class CathSql(SqlBuilder):

    def __init__(self, desc_file, sql_db, table_name, column_definition):
        super(CathSql, self).__init__(sql_db, table_name, column_definition)

        self.file = desc_file
        self.records = []

    def parse_desc(self):
        record_items = {}
        prefixes = ["DOMAIN", "CATHCODE", "CLASS", "ARCH", "TOPOL", "HOMOL"]
        reg = re.compile(r"  +")

        for line in self.read_text_generator(self.file):
            for prefix in prefixes:
                if line.startswith(prefix):
                    record_items[prefix] = reg.split(line)[1].rstrip()
            if line.startswith("//\n"):
                self.records.append(tuple(record_items.values()))

    def run(self):
        logger.info(f"Start to create SQL table: {table_name} in SQL file {sql_db}")

        logger.info(f"Step 1 -> Create SQL table: {table_name}")
        cursor = self.create_cath_sql()

        logger.info("Step 2 -> Format records for SQL batch inserting")
        self.parse_desc()

        logger.info(f"Step 3 -> Insert data in bulk to SQL table: {table_name}")
        self.insert_batch(cursor, self.records, 6)

        logger.info(f'Step 4 -> Creating index for SQL table: {table_name}')
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
    cath_builder = CathSql(desc_file, sql_db, table_name, column_definition)
    cath_builder.run()
