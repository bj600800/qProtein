"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/09/13

# Description: Parse blast.
# ------------------------------------------------------------------------------
"""
from qprotein.preprocessing.blast import BlastQuery
from qprotein.preprocessing.sqlite3_builder import SqlBuilder
from qprotein.preprocessing.sqlite3_searcher import SqlSearch
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)

class SourceMap(SqlBuilder, SqlSearch):
	def __init__(self, sql_db, pdb_blast_out, uniprot_blast_out):
		super().__init__()
		self.pdb_blast_out = pdb_blast_out
		self.uniprot_blast_out = uniprot_blast_out
		self.sql_db = sql_db
		self.table_name = "results_summary"
		self.column_definition = [
			('query_name', 'TEXT'), ('source', 'TEXT'), ('accession_id', 'TEXT'),
			('subject_start', 'INT'), ('subject_end', 'INT')]
		
	def format_records(self, hit, source):
		record_items = {column_name[0]: '' for column_name in self.column_definition}
		record_items['query_name'] = hit[0]
		record_items['source'] = source
		if source == "AFDB":
			accession = hit[1].split("|")[1]
		else:
			accession = hit[1]
		record_items['accession_id'] = accession
		record_items['subject_start'] = hit[2]
		record_items['subject_end'] = hit[3]
		record = [tuple(record_items.values())]
		return record
	
	def parse_blast_out(self, cursor, columns):
		total_record_count = 0
		for i, blast_out in enumerate([self.pdb_blast_out, self.uniprot_blast_out]):
			records = []
			record_count = 0
			source = "PDB" if i == 0 else "AFDB"
			if blast_out:
				for hit in blast_out:
					record = self.format_records(hit.split("\t"), source)
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
		self.parse_blast_out(cursor=cursor, columns=columns)
		self.create_index(cursor=cursor, table_name=self.table_name, columns=columns)

if __name__ == '__main__':
	sql_db = r"D:\subject\active\1-qProtein\data\enzymes\GH8\1_preprocessing\qprotein_results.db"
	db_dir = r"D:\subject\active\1-qProtein\database"
	blast = BlastQuery(sql_db=sql_db, db_dir=db_dir)
	pdb_blast_out, uniprot_blast_out = blast.run()
	mapper = SourceMap(sql_db=sql_db, pdb_blast_out=pdb_blast_out, uniprot_blast_out=uniprot_blast_out)
	mapper.run()
	
