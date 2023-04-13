"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: parse uniprot dat file to sqlite3
# ------------------------------------------------------------------------------
"""


from qprotein.database.sqlite3_builder import SqlBuilder
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


class UniprotSql(SqlBuilder):
    def __init__(self, dat_file, sql_db, table_name, column_definition):
        super(UniprotSql, self).__init__(sql_db, table_name, column_definition)

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

    def parse_dat(self, cursor):
        for line in self.read_gz_generator(self.dat_file):
            if self.record_counter == 500000:
                self.insert_many(cursor, self.records, 8)
                logger.info(f"Insert 500000 records into {self.table_name}, total records:" + str(self.total_records))

                self.records = []
                self.record_counter = 0

            else:
                self.format_records(line)

        # Insert the remaining records
        if self.record_counter > 0:
            logger.info(f"Insert {self.record_counter} records into {self.table_name}")
            self.insert_many(cursor, self.records, 8)

        logger.info(f"Total records:" + str(self.total_records))

    def run(self):
        columns = [i[0] for i in self.column_definition]
        logger.info(f"Start to create SQL table: {self.table_name} in SQL file {self.sql_path}")
        logger.info(f"Create SQL table: {self.table_name}")
        cursor = self.create_table(self.table_name, self.sql_path, self.column_definition)

        logger.info(f"Parse uniprot dat file and insert records into {self.table_name}")
        self.parse_dat(cursor)

        logger.info(f'Create index for SQL table: {self.table_name}')
        self.create_index(cursor=cursor, columns=columns)

        logger.info(f'Successfully built SQL database for {self.table_name}')


if __name__ == '__main__':
    sprot_dat = r'G:\DB\uniprot_sprot.dat.gz'
    sql_db = r'C:\Users\bj600\Desktop\qprotein_db.db'
    table_name = 'sprot_dat'
    column_definition = {"accession": "TEXT", "gene_name": "TEXT",
                         "ec_number": "TEXT", "go_component": "TEXT",
                         "go_process": "TEXT", "go_function": "TEXT",
                         "interpro": "TEXT", "pfam": "TEXT"
                         }

    uniprot_builder = UniprotSql(sprot_dat, sql_db, table_name, column_definition)
    uniprot_builder.run()
