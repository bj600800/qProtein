"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/09/15

# Description: Update structure path to sql results table
# ------------------------------------------------------------------------------
"""

from qprotein.preprocessing.sqlite3_builder import SqlBuilder
from qprotein.preprocessing.sqlite3_searcher import SqlSearch
from qprotein.utilities import logger
logger = logger.setup_log(name=__name__)

class UpdateSQL(SqlBuilder, SqlSearch):
	def __init__(self, sql_db, struct_info):
		super().__init__()
		self.sql_db = sql_db
		self.struct_info = struct_info
		self.table_name = "results_summary"
		self.column_definition = [('query_name', 'TEXT'), ('struct_path', 'TEXT')]
	
	def format_records(self, info):
		record_items = {column_name[0]: '' for column_name in self.column_definition}
		record_items['query_name'] = info[0]
		record_items['struct_path'] = info[1]
		record = [tuple(record_items.values())]
		return record
	
	def parse_struct_info(self, cursor, columns):
		total_record_count = 0
		for info in self.struct_info:
			records = []
			record_count = 0
			record = self.format_records(info)
			records.extend(record)
			record_count += 1
			total_record_count += 1
			if record_count == 500000:
				self.update_many_columns(
						cursor=cursor,
						columns=columns,
						records=records,
						table_name=self.table_name
						)
				records = []
				record_count = 0
				
			# Insert the remaining records
			if record_count > 0:
				self.update_many_columns(
						cursor=cursor,
						columns=columns,
						records=records,
						table_name=self.table_name
						)
		logger.info(f"Update records: {total_record_count}")
	
	def run(self):
		columns = [i[0] for i in self.column_definition]
		cursor = self.connect_sql(self.sql_db)
		self.add_column(cursor=cursor, table_name="results_summary", column_definition=self.column_definition)
		self.parse_struct_info(cursor=cursor, columns=columns)
		self.create_index(cursor=cursor, table_name=self.table_name, columns=columns)

import os
from qprotein.mapping.PDB import PDBMapper
from qprotein.mapping.AFDB import AFDBMapper

if __name__ == '__main__':
	structure_folder = r"D:\subject\active\1-qProtein\data\enzymes\GH8\2_StrucMapping"
	sql_db = r"D:\subject\active\1-qProtein\data\enzymes\GH8\1_preprocessing\qprotein_results.db"
	info = []
	PDB = PDBMapper(sql_db=sql_db, folder=structure_folder)
	pdb_info = PDB.run()
	AFDB = AFDBMapper(sql_db=sql_db, folder=structure_folder)
	afdb_info = AFDB.run()
	for d in [os.path.join(structure_folder, "positive"), os.path.join(structure_folder, "negative")]:
		for file in os.listdir(d):
			struc_path = os.path.join(d, file)
			struc_name = os.path.splitext(file)[0]
			info.append([struc_name, struc_path])
		update = UpdateSQL(sql_db=sql_db, struct_info=info)
		update.run()
	