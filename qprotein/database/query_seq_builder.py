"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/03/13

# Description: create sql database for users
# ------------------------------------------------------------------------------
"""
from qprotein.utilities import logger
from qprotein.database.sqlite3_builder import SqlBuilder
from Bio import SeqIO

logger = logger.setup_log(name=__name__)


class QuerySql(SqlBuilder):
    def __init__(self, sql_path, fasta_file):
        super(QuerySql, self).__init__()
        self.sql_db = sql_path
        self.table_name = "results_summary"
        self.column_definition = [("query_name", "TEXT"), ("sequence", "TEXT")]
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
                self.insert_many_columns(cursor=cursor,
                                         values_num=2,
                                         columns=columns,
                                         records=self.records,
                                         table_name=self.table_name
                                         )
                logger.info(f"Insert 500000 records. Total records: {str(self.total_records)}")
                self.records = []
                self.record_counter = 0
                self.format_records(seq_record)
            else:
                self.format_records(seq_record)

        # Insert the remaining records
        if self.record_counter > 0:
            self.insert_many_columns(cursor=cursor,
                                     values_num=2,
                                     columns=columns,
                                     records=self.records,
                                     table_name=self.table_name
                                     )
            logger.info(f"Insert {self.record_counter} records")
        logger.info(f"Total records:" + str(self.total_records))

    def run(self):
        logger.info(f"Parse {self.table_name} sequences and insert records")
        columns = [i[0] for i in self.column_definition]
        sql = f"CREATE TABLE {self.table_name} (query_name TEXT PRIMARY KEY, sequence TEXT)"
        cursor = self.create_table(table_name=self.table_name, sql_db=self.sql_db, sql=sql)
        self.parse_fasta(cursor=cursor, columns=columns)
        self.create_index(cursor=cursor, table_name=self.table_name, columns=columns)
        logger.info(f'Successfully built SQL database for task: {self.table_name}')

