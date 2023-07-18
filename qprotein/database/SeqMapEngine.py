"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: search the sequences with existing structures for query sequences from PDB & Alphafold databases
# ------------------------------------------------------------------------------
"""

from qprotein.database.sqlite3_builder import SqlBuilder
from qprotein.database.sqlite3_searcher import SqlSearch
from qprotein.database.query_seq_builder import QuerySql
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


def create_user_table(summary_sql_path, fasta_file):
	creater = QuerySql(sql_path=summary_sql_path, fasta_file=fasta_file)
	creater.run()
	
class SprotDmnd(SqlBuilder, SqlSearch):
	"""
    find the same structure of query sequence from sprot sequence database
    """
	
	def __init__(self, dmnd_output, summary_sql_path):
		super().__init__()
		self.dmnd_output = dmnd_output
		self.summary_sql_path = summary_sql_path
		self.idx_name_prefix = 'sprot'
		self.table_name = "results_summary"
		self.column_definition = [
			('query_name', 'TEXT'),
			('query_length', 'TEXT'),
			('sprot_acc', 'TEXT'),
			('sprot_subject_start', 'INT'),
			('sprot_subject_end', 'INT')
			]
		self.record_items = {column_name[0]: '' for column_name in self.column_definition}
		self.records = []
		self.record_counter = 0
		self.total_records = 0
	
	def parse_items_from_dmnd(self, line):
		self.record_items['query_name'] = line[0]
		self.record_items['query_length'] = line[5]
		self.record_items['sprot_acc'] = line[1].split("|")[1]
		self.record_items['sprot_subject_start'] = line[7]
		self.record_items['sprot_subject_end'] = line[8]
	
	def format_dmnd(self, line):
		# prepare
		self.parse_items_from_dmnd(line)
		self.records.append(tuple(self.record_items.values()))
		self.record_counter += 1
		self.total_records += 1
	
	def operate_dmnd(self, cursor, columns):
		for line in self.read_text_generator(filename=self.dmnd_output, header=False):
			line = line.split('\t')
			if self.record_counter == 500000:
				self.insert_update_columns(
						cursor=cursor, table_name="results_summary", columns=columns,
						records=self.records
						)
				logger.info(f"Insert 500000 records. Total records: {str(self.total_records)}")
				self.records = []
				self.record_counter = 0
			else:
				self.format_dmnd(line)
		
		if self.record_counter > 0:
			logger.info(f"Insert {self.record_counter} records")
			self.insert_update_columns(
					cursor=cursor, table_name="results_summary", columns=columns,
					records=self.records
					)
		
		logger.info(f"Total records:" + str(self.total_records))
	
	def run(self):
		logger.info(f"Parse {self.idx_name_prefix} diamond outputs and insert records")
		columns = [i[0] for i in self.column_definition]
		cursor = self.connect_sql(sql_db=self.summary_sql_path)
		self.add_column(cursor=cursor, table_name="results_summary", column_definition=self.column_definition)
		self.operate_dmnd(cursor=cursor, columns=columns)
		self.create_index(cursor=cursor, columns=columns, table_name=self.table_name)


class TremblDmnd(SprotDmnd):
	"""
    find the same structure of the query sequence from trembl sequence database
    """
	
	def __init__(self, dmnd_output, summary_sql_path):
		super().__init__(
				dmnd_output=dmnd_output,
				summary_sql_path=summary_sql_path
				)
		self.column_definition = [
			('query_name', 'TEXT'),
			('query_length', 'TEXT'),
			('trembl_acc', 'TEXT'),
			('trembl_subject_start', 'INT'),
			('trembl_subject_end', 'INT')
			]
		self.idx_name_prefix = 'trembl'
		self.record_items = {column_name[0]: '' for column_name in self.column_definition}
		self.records = []
		self.record_counter = 0
		self.total_records = 0
	
	def parse_items_from_dmnd(self, line):
		self.record_items['query_name'] = line[0]
		self.record_items['query_length'] = line[5]
		self.record_items['trembl_acc'] = line[1].split("|")[1]
		self.record_items['trembl_subject_start'] = line[7]
		self.record_items['trembl_subject_end'] = line[8]


import os
if __name__ == '__main__':
	work_path = r"D:\subject\active\1-qProtein\data\enzymes\cazy_endo-1_4-beta-xylanase"
	summary_sql_path = os.path.join(work_path, "qprotein_results.db")
	dmnd_output = os.path.join(work_path, "17_length_rmdup_temperature_cazy_domain_alphafolddb_trembl_filtered.out")
	fasta_file = os.path.join(work_path, "9_length_rmdup_temperature_cazy_domain.fasta")
	create_user_table(summary_sql_path=summary_sql_path, fasta_file=fasta_file)
	create_sprot = TremblDmnd(dmnd_output=dmnd_output, summary_sql_path=summary_sql_path)
	create_sprot.run()
	
