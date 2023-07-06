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
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


class CreateSql(SqlBuilder):
	def __init__(self, summary_sql_path):
		super(SqlBuilder, self).__init__()
		self.summary_sql_path = summary_sql_path
		self.create()
	
	def create(self):
		sql = f"CREATE TABLE 'results_summary' (query_name TEXT PRIMARY KEY)"
		self.create_table(table_name='results_summary', sql_db=self.summary_sql_path, sql=sql)


class SprotDmnd(SqlBuilder, SqlSearch):
	"""
    find the same structure of query sequence from sprot sequence database
    """
	
	def __init__(self, dmnd_output, summary_sql_path):
		super().__init__()
		self.sql_db = summary_sql_path,
		self.table_name = "results_summary"
		self.column_definition = [
			('query_name', 'TEXT'),
			('sprot_acc', 'TEXT'),
			('sprot_subject_start', 'INT'),
			('sprot_subject_end', 'INT')
			]
		self.dmnd_output = dmnd_output
		self.summary_sql_path = summary_sql_path
		self.idx_name_prefix = 'sprot'
		self.record_items = {column_name[0]: '' for column_name in self.column_definition}
		self.records = []
		self.record_counter = 0
		self.total_records = 0
	
	def parse_items_from_dmnd(self, line):
		self.record_items['query_name'] = line[0]
		self.record_items['sprot_subject_start'] = line[6]
		self.record_items['sprot_subject_end'] = line[7]
	
	def format_dmnd(self, line):
		# prepare
		self.parse_items_from_dmnd(line)
		self.records.append(tuple(self.record_items.values()))
		self.record_counter += 1
		self.total_records += 1
	
	def parse_dmnd(self, cursor, columns):
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
		self.parse_dmnd(cursor=cursor, columns=columns)
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
			('trembl_acc', 'TEXT'),
			('trembl_subject_start', 'INT'),
			('trembl_subject_end', 'INT'),
			]
		self.idx_name_prefix = 'trembl'
		self.record_items = {column_name[0]: '' for column_name in self.column_definition}
		self.records = []
		self.record_counter = 0
		self.total_records = 0
	
	def parse_items_from_dmnd(self, line):
		self.record_items['query_name'] = line[0]
		self.record_items['trembl_acc'] = line[1].split('|')[1]
		self.record_items['trembl_subject_start'] = line[6]
		self.record_items['trembl_subject_end'] = line[7]


class QueryAnalysis(SqlBuilder, SqlSearch):
	"""
    analyze the query sequences and insert the results to sql
    """
	
	def __init__(self, task_name, summary_sql_path, fasta_sql_path):
		super().__init__()
		self.sql_db = summary_sql_path
		self.table_name = "results_summary"
		self.column_definition = [('query_name', 'TEXT'), ('seq_length', 'INT')]
		self.summary_sql_path = summary_sql_path
		self.fasta_sql_path = fasta_sql_path
		self.record_items = {column_name[0]: '' for column_name in self.column_definition}
		self.task_name = task_name
		self.total_records = 0
		
	def get_query_name(self):
		summary_cursor = self.connect_sql(sql_db=self.summary_sql_path)
		sql_cmd = "SELECT query_name FROM results_summary WHERE sprot_acc != '' or trembl_acc != ''"
		query_name = [i[0] for i in SqlSearch.fetch_results(cursor=summary_cursor, sql_cmd=sql_cmd)]
		return query_name
	
	def get_query_length(self):
		fasta_cursor = self.connect_sql(sql_db=self.fasta_sql_path)
		query_name = self.get_query_name()
		batch_size = 100000
		query_length = []
		total_query = len(query_name)
		for i in range(0, total_query, batch_size):
			batch_acc = query_name[i:i + batch_size]
			sql_cmd = "SELECT query_name, length(sequence) FROM {} WHERE query_name IN ({})" \
				.format(self.task_name, ', '.join(['?'] * len(batch_acc)))
			return_dat = SqlSearch.fetch_many_results(cursor=fasta_cursor, sql_cmd=sql_cmd, search_targets=batch_acc)
			query_length.extend(return_dat)
		return query_length
	
	def parse_items_from_db(self, dat):
		self.record_items['query_name'] = dat[0]
		self.record_items['seq_length'] = dat[1]
	
	def format_records(self, dat):
		self.parse_items_from_db(dat)
		self.records.append(tuple(self.record_items.values()))
		self.record_counter += 1
		self.total_records += 1
	
	def parse_uniprot(self, cursor, columns):
		query_length = self.get_query_length()
		for i in query_length:
			if self.record_counter == 500000:
				self.update_many_columns(
						cursor=cursor, table_name="results_summary",
						columns=columns, records=self.records
						)
				logger.info(f"Insert 500000 records. Total records: {str(self.total_records)}")
				
				self.records = []
				self.record_counter = 0
			else:
				self.format_records(i)
		if self.record_counter > 0:
			self.update_many_columns(cursor=cursor, table_name="results_summary", columns=columns, records=self.records)
			logger.info(f"Update {self.record_counter} records")
		logger.info(f"Total records:" + str(self.total_records))
	
	def run(self):
		logger.info(f"Add query_length for each sequence")
		columns = [i[0] for i in self.column_definition]
		cursor = self.connect_sql(sql_db=self.summary_sql_path)
		self.add_column(
				cursor=cursor,
				table_name="results_summary",
				column_definition=self.column_definition[1:],
				specific_column='query_name'
				)
		self.parse_uniprot(cursor=cursor, columns=columns)
		self.create_index(cursor=cursor, columns=columns, table_name=self.table_name)
		