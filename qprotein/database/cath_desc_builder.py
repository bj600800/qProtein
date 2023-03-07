"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: Build cath description sqlite file, source file read in by a generator from SqlBuilder
# ------------------------------------------------------------------------------
"""

import re

from qprotein.utilities import logger
from qprotein.database.sqlite3_builder import SqlBuilder

logger = logger.setup_log(name=__name__)


class CathSql(SqlBuilder):

    def __init__(self, desc_file, sql_db, table_name, column_definition):
        super(CathSql, self).__init__(sql_db, table_name, column_definition)

        self.file = desc_file

        self.record_items = {}
        self.records = []
        self.record_counter = 0
        self.total_records = 0

    def format_records(self, line):
        prefixes = ["DOMAIN", "CATHCODE", "CLASS", "ARCH", "TOPOL", "HOMOL"]
        reg = re.compile(r"  +")
        for prefix in prefixes:
            if line.startswith(prefix):
                self.record_items[prefix] = reg.split(line)[1].rstrip()

        if line.startswith("//\n"):
            self.records.append(tuple(self.record_items.values()))
            self.record_counter += 1
            self.total_records += 1

    def parse_desc(self, cursor):
        for line in self.read_text_generator(self.file, header=False):
            if self.record_counter == 500000:
                logger.info(f"Insert 500000 records into {self.table_name}, total records:" + str(self.total_records))
                self.insert_many(cursor, self.records, 6)
                self.records = []
                self.record_counter = 0

            else:
                self.format_records(line)

        if self.record_counter > 0:
            logger.info(f"Insert {self.record_counter} records into {self.table_name}")
            self.insert_many(cursor, self.records, 6)

    def run(self):
        logger.info(f"Start to create SQL table: {self.table_name} in SQL file {self.sql_db}")

        logger.info(f"Create SQL table: {self.table_name}")
        cursor = self.create_table(self.table_name, self.sql_db, self.column_definition)

        logger.info(f"Parse file and insert records into {self.table_name}")
        self.parse_desc(cursor)

        logger.info(f'Creating index for SQL table: {self.table_name}')
        self.create_index_all(cursor)

        logger.info(f'Successfully built SQL database for {self.table_name}')




if __name__ == '__main__':
    desc_file = r'G:\DB\cath\cath-classification-data_txt\cath-domain-description-file.txt'
    sql_db = r'C:\Users\bj600\Desktop\qprotein_db.db'
    table_name = 'cath_desc'
    column_definition = {"domain_name": "TEXT", "cath_code": "TEXT",
                         "cath_class": "TEXT", "cath_architecture": "TEXT",
                         "cath_topology": "TEXT", "cath_homologous": "TEXT"
                         }
    cath_builder = CathSql(desc_file, sql_db, table_name, column_definition)
    cath_builder.run()
