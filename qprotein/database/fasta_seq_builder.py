"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/13

# Description: 
# ------------------------------------------------------------------------------
"""

from qprotein.utilities import logger
from qprotein.database.sqlite3_builder import SqlBuilder
from Bio import SeqIO

logger = logger.setup_log(name=__name__)


class FastaSql(SqlBuilder):
    def __init__(self, fasta_file, sql_db, table_name, column_definition):
        super(FastaSql, self).__init__(sql_db, table_name, column_definition)

        self.fasta_file = fasta_file

        self.record_items = {column_name[0]: '' for column_name in self.column_definition}
        self.records = []
        self.record_counter = 0
        self.total_records = 0

    def parse_items(self, seq_record):
        self.record_items['query_name'] = str(seq_record.id)
        self.record_items['sequence'] = str(seq_record.seq)

    def format_records(self, seq_record):
        self.parse_items(seq_record)
        self.records.append(tuple(self.record_items.values()))
        self.record_counter += 1
        self.total_records += 1

    def parse_fasta(self, cursor):
        for seq_record in SeqIO.parse(self.fasta_file, "fasta"):
            if self.record_counter == 500000:
                self.insert_many(cursor, self.records, 2)
                logger.info(f"Insert 500000 records into {self.table_name}, total records:" + str(self.total_records))

                self.records = []
                self.record_counter = 0

            else:
                self.format_records(seq_record)

            # Insert the remaining records
        if self.record_counter > 0:
            logger.info(f"Insert {self.record_counter} records into {self.table_name}")
            self.insert_many(cursor, self.records, 2)

        logger.info(f"Total records:" + str(self.total_records))

    def run(self):
        columns = [i[0] for i in self.column_definition]

        logger.info(f"Start to create SQL table: {self.table_name} in SQL file {self.sql_db}")

        logger.info(f"Create SQL table: {self.table_name}")
        sql = f"CREATE TABLE {self.table_name} (query_name TEXT PRIMARY KEY, sequence TEXT)"
        cursor = self.create_table(self.table_name, self.sql_db, sql)

        logger.info(f"Parse file and insert records into {self.table_name}")
        self.parse_fasta(cursor)

        logger.info(f'Creating index for SQL table: {self.table_name}')
        self.create_index(cursor, columns)

        logger.info(f'Successfully built SQL database for {self.table_name}')


if __name__ == '__main__':
    fasta_file = r'D:\subject\active\1-qProtein\data\manure\manure.fasta'
    sql_db = r'D:\subject\active\1-qProtein\data\manure\qprotein_db.db'
    table_name = 'query_seq'
    column_definition = [("query_name", "TEXT"), ("sequence", "TEXT")]
    fasta_db = FastaSql(fasta_file, sql_db, table_name, column_definition)
    fasta_db.run()
