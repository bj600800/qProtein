"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/03/13

# Description: Build sql table for the input fasta.

# ------------------------------------------------------------------------------
"""
import os.path
import re
from Bio import SeqIO

from qprotein.preprocessing.sqlite3_builder import SqlBuilder
from qprotein.utilities import logger
logger = logger.setup_log(name=__name__)


class BuildSql(SqlBuilder):
    def __init__(self, sql_path, positive_file, negative_file):
        super(BuildSql, self).__init__()
        self.sql_db = sql_path
        self.table_name = "results_summary"
        self.column_definition = [("query_name", "TEXT"), ("label", "TEXT"), ("sequence", "TEXT")]
        self.positive_file = positive_file
        self.negative_file = negative_file
        self.run()
        
    def format_records(self, seq_record, label):
        record_items = {column_name[0]: '' for column_name in self.column_definition}
        record_items['query_name'] = str(seq_record.id)
        record_items['label'] = label
        record_items['sequence'] = str(seq_record.seq)
        record = [tuple(record_items.values())]
        return record

    def parse_fasta(self, cursor, columns):
        total_record_count = 0
        for i, fasta_file in enumerate([self.positive_file, self.negative_file]):
            records = []
            record_count = 0
            label = "positive" if i == 0 else "negative"
            for seq_record in SeqIO.parse(fasta_file, "fasta"):
                special_chars = ["[", "]", "/", ",", "#", "&", "\\", "(", ")", "'", "."," "]
                replace_char = "_"
                for char in special_chars:
                    seq_record.id = seq_record.id.replace(char, replace_char)
                seq_record.id = re.sub("_+", "_", seq_record.id)
                record = self.format_records(seq_record, label)
                records.extend(record)
                record_count += 1
                total_record_count += 1

                if record_count == 500000:
                    self.insert_many_columns(cursor=cursor,
                                             values_num=3,
                                             columns=columns,
                                             records=records,
                                             table_name=self.table_name
                                             )
                    records = []
                    record_count = 0

            # Insert the remaining records
            if record_count > 0:
                self.insert_many_columns(cursor=cursor,
                                         values_num=3,
                                         columns=columns,
                                         records=records,
                                         table_name=self.table_name
                                         )
        logger.info(f"Total records: {total_record_count}")

    def run(self):
        columns = [i[0] for i in self.column_definition]
        sql = f"CREATE TABLE {self.table_name} (query_name TEXT PRIMARY KEY, label TEXT, sequence TEXT)"
        cursor = self.create_table(table_name=self.table_name, sql_db=self.sql_db, sql=sql)
        self.parse_fasta(cursor=cursor, columns=columns)
        self.create_index(cursor=cursor, table_name=self.table_name, columns=columns)


if __name__ == '__main__':
    wd = input("Data absolute directory: ")  # D:\subject\active\1-qProtein\code\test
    sql_path = os.path.join(wd, "qprotein_results.db")
    positive_filename = input("Positive fasta filename: ")  # positive.fasta
    negative_filename = input("Negative fasta filename: ")  # negative.fasta
    positive_fasta = os.path.join(wd, positive_filename)
    negative_fasta = os.path.join(wd, negative_filename)
    test = BuildSql(sql_path, positive_fasta, negative_fasta)
