"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: create sql database for CATH description
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

        self.record_items = {column_name[0]: '' for column_name in self.column_definition}
        self.records = []
        self.record_counter = 0
        self.total_records = 0

    def format_records(self, line):
        reg = re.compile(r"  +")
        if line.startswith("DOMAIN"):
            self.record_items["domain_name"] = reg.split(line)[1].rstrip()
        elif line.startswith("CATHCODE"):
            self.record_items["cath_code"] = reg.split(line)[1].rstrip()
        elif line.startswith("CLASS"):
            self.record_items["cath_class"] = reg.split(line)[1].rstrip()
        elif line.startswith("ARCH"):
            self.record_items["cath_architecture"] = reg.split(line)[1].rstrip()
        elif line.startswith("TOPOL"):
            self.record_items["cath_topology"] = reg.split(line)[1].rstrip()
        elif line.startswith("HOMOL"):
            self.record_items["cath_homologous"] = reg.split(line)[1].rstrip()
        elif line.startswith("//\n"):
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
        columns = [i[0] for i in self.column_definition]
        logger.info(f"Start to create SQL table: {self.table_name} in SQL file {self.sql_path}")

        logger.info(f"Create SQL table: {self.table_name}")
        sql = f"CREATE TABLE {self.table_name} ({', '.join([i[0]+' '+i[1] for i in self.column_definition])})"
        cursor = self.create_table(self.table_name, self.sql_path, sql)

        logger.info(f"Parse file and insert records into {self.table_name}")
        self.parse_desc(cursor)

        logger.info(f'Creating index for SQL table: {self.table_name}')
        self.create_index(cursor, columns)

        logger.info(f'Successfully built SQL database for {self.table_name}')


if __name__ == '__main__':
    desc_file = r'G:\DB\cath\cath-classification-data_txt\cath-domain-description-file.txt'
    sql_db = r'D:\subject\active\1-qProtein\data\tibet\qprotein_db.db'
    table_name = 'cath_desc'
    column_definition = [("domain_name", "TEXT"), ("cath_code", "TEXT"),
                         ("cath_class", "TEXT"), ("cath_architecture", "TEXT"),
                         ("cath_topology", "TEXT"), ("cath_homologous", "TEXT")]

    cath_builder = CathSql(desc_file, sql_db, table_name, column_definition)
    cath_builder.run()
