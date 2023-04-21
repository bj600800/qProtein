"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: parse uniprot dat file with sqlite3
# ------------------------------------------------------------------------------
"""

import os
from tqdm import tqdm
from qprotein.database.sqlite3_builder import SqlBuilder
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


class UniprotSql(SqlBuilder):
    def __init__(self, dat_file, sql_path, table_name):
        super(UniprotSql, self).__init__(sql_db=sql_path, table_name=table_name,
                                         column_definition=[("accession", "TEXT"), ("gene_name", "TEXT"),
                                                            ("ec_number", "TEXT"), ("go_component", "TEXT"),
                                                            ("go_process", "TEXT"), ("go_function", "TEXT"),
                                                            ("interpro", "TEXT"), ("pfam", "TEXT")]
                                         )

        self.dat_file = dat_file

        self.record_items = {'accession': '', 'geneName': '', 'ec_number': [],
                             'go_cellular': [], 'go_process': [], 'go_function': [],
                             'interpro': [], 'pfam': []
                             }

        self.records = []
        self.record_counter = 0
        self.total_records = 0

    def parse_accession(self, line):
        return line.split()[1].replace(';', '')

    def parse_gene_name(self, line):
        return line.split("Full=")[1].split("{")[0].rstrip().lower()

    def parse_ec_code(self, line):
        ec_code = line.strip().split()[1].replace(";", "").replace("EC=", "").lstrip()
        return ec_code

    def parse_DR(self, line):
        return '|'.join(line.rstrip().split('; ')[1:3]).rstrip('.')

    def parse_go(self, line):
        if "; C:" in line:
            go = self.parse_DR(line)
            self.record_items['go_cellular'].append(go)

        elif "; F:" in line:
            go = self.parse_DR(line)
            self.record_items['go_function'].append(go)

        elif "; P:" in line:
            # biological process
            go = self.parse_DR(line)
            self.record_items['go_process'].append(go)

    def parse_end_mark(self):
        self.records.append(tuple(self.format_items(self.record_items).values()))

        # init items
        self.record_items = {'accession': '', 'geneName': '', 'ec_number': [],
                             'go_cellular': [], 'go_process': [], 'go_function': [],
                             'interpro': [], 'pfam': []
                             }

    def format_items(self, record_items):
        format_output = {}
        for k, v in record_items.items():
            if type(v) == list:
                v = ';'.join(v)
            format_output[k] = v
        return format_output

    def format_records(self, line):
        if line.startswith("AC"):
            self.record_items['accession'] = self.parse_accession(line)

        elif line.startswith("DE   RecName"):
            self.record_items['geneName'] = self.parse_gene_name(line)

        elif line.startswith("DE") and "EC=" in line:
            ec_code = self.parse_ec_code(line)
            self.record_items['ec_number'].append(ec_code)

        elif line.startswith("DR   InterPro"):
            ipr = self.parse_DR(line)
            self.record_items['interpro'].append(ipr)

        elif line.startswith("DR   Pfam"):
            pf = self.parse_DR(line)
            self.record_items['pfam'].append(pf)

        elif line.startswith("DR   GO"):
            self.parse_go(line)

        elif line.startswith("//\n"):
            self.parse_end_mark()
            self.record_counter += 1
            self.total_records += 1

    def parse_dat(self, cursor, columns):
        for line in self.read_gz_generator(self.dat_file):
            if self.record_counter == 500000:
                self.insert_update_columns(cursor=cursor, table_name=self.table_name, columns=columns,
                                           records=self.records
                                           )
                logger.info(f"Insert 500000 records. Total records: {str(self.total_records)}")

                self.records = []
                self.record_counter = 0

            else:
                self.format_records(line)

        # Insert the remaining records
        if self.record_counter > 0:
            self.insert_update_columns(cursor=cursor, table_name=self.table_name, columns=columns,
                                       records=self.records
                                       )
            logger.info(f"Insert {self.record_counter} records")
        logger.info(f"Total records:" + str(self.total_records))

    def run(self):
        logger.info(f"Parse {self.table_name} and insert records")
        columns = [i[0] for i in self.column_definition]

        sql = f"CREATE TABLE {self.table_name} (accession TEXT PRIMARY KEY)"
        cursor = self.create_table(self.table_name, self.sql_path, sql)
        self.add_column(cursor=cursor, table_name=self.table_name, column_definition=self.column_definition)
        self.parse_dat(cursor=cursor, columns=columns)
        self.create_index(cursor=cursor, columns=columns)
        logger.info(f'Successfully built SQL database for task: {self.table_name}')
