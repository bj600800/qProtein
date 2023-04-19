"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/13

# Description: create sql database for fasta sequences
# ------------------------------------------------------------------------------
"""
import os.path

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

    def parse_fasta(self, cursor, columns):
        for seq_record in SeqIO.parse(self.fasta_file, "fasta"):
            if self.record_counter == 500000:
                self.insert_many_columns(cursor=cursor, values_num=2,
                                         columns=columns, records=self.records
                                         )
                logger.info(f"Insert 500000 records. Total records:" + str(self.total_records))
                self.records = []
                self.record_counter = 0
                self.format_records(seq_record)
            else:
                self.format_records(seq_record)

            # Insert the remaining records
        if self.record_counter > 0:
            self.insert_many_columns(cursor=cursor, values_num=2,
                                     columns=columns, records=self.records
                                     )
            logger.info(f"Update {self.record_counter} records")

        logger.info(f"Total records:" + str(self.total_records))

    def run(self):
        columns = [i[0] for i in self.column_definition]

        logger.info(f"Start to create SQL table: {self.table_name} in SQL file {self.sql_path}")

        logger.info(f"Create SQL table: {self.table_name}")
        sql = f"CREATE TABLE {self.table_name} (query_name TEXT PRIMARY KEY, sequence TEXT)"
        cursor = self.create_table(self.table_name, self.sql_path, sql)

        logger.info(f"Parse file and insert records into {self.table_name}")
        self.parse_fasta(cursor, columns)

        logger.info(f'Creating index for SQL table: {self.table_name}')
        self.create_index(cursor, columns)

        logger.info(f'Successfully built SQL database for {self.table_name}')


if __name__ == '__main__':
    task_name = 'tibet'
    root_dir = r'D:\subject\active\1-qProtein\data'
    work_dir = os.path.join(root_dir, task_name)
    fasta_file = os.path.join(work_dir, 'tibet.fasta')
    sql_db = os.path.join(root_dir, 'qprotein_db.db')
    table_name = os.path.split(fasta_file)[1].split('.')[0]
    column_definition = [("query_name", "TEXT"), ("sequence", "TEXT")]
    fasta_db = FastaSql(fasta_file, sql_db, table_name, column_definition)
    fasta_db.run()