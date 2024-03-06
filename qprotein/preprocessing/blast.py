"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/09/09

# Description: Run blast and insert results

# ------------------------------------------------------------------------------
"""
import os
import tempfile
from diamond4py import Diamond, OutFormat

from qprotein.preprocessing.sqlite3_searcher import SqlSearch
from qprotein.preprocessing.sqlite3_builder import SqlBuilder
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


def seq2file(query_fasta):
	with tempfile.NamedTemporaryFile(mode='w', delete=False) as query_file:
		for seq in query_fasta:
			query_file.write(">" + seq[0] + "\n")
			query_file.write(seq[1] + "\n")
		query_file = query_file.name
	return query_file

def read_out(out_file):
	with open(out_file, "r") as f:
		content = f.readlines()
	out = [line.rstrip().split("\t") for line in content]
	return out
	

class RunBlast(SqlSearch, SqlBuilder):
	def __init__(self, wd, sql_db, db_dir):
		self.wd = wd
		self.sql_db = sql_db
		self.db_dir = db_dir
		self.run()
		
	def get_seq(self, sql_cmd, search_targets=None):
		cursor = self.connect_sql(self.sql_db)
		if search_targets:
			query_fasta = self.fetch_many_results(cursor, sql_cmd, search_targets)
		else:
			query_fasta = self.fetch_results(cursor, sql_cmd)
		return query_fasta
	
	def blastp(self, wd, query_fasta, db):
		query_file = seq2file(query_fasta)
		diamond = Diamond(
				database=db,
				n_threads=4
				)
		OutFormat.BLAST_TABULAR.with_extra_option(
				"qseqid","sseqid","sstart", "send"
				)
		out_file = os.path.join(wd, os.path.basename(db).split(".")[0]+'.out')
		diamond.blastp(
				query=query_file,
				out=out_file,
				outfmt=OutFormat.BLAST_TABULAR,
				sensitivity=4,
				max_target_seqs=1,
				evalue=1e-10,
				id=100,
				query_cover=100)
		os.remove(query_file)
		return out_file
	

	def run_uniprot(self, query_id, dmnd):
		cmd = "SELECT query_name, sequence FROM results_summary WHERE query_name IN ({})" \
			.format(', '.join(['?'] * len(query_id)))
		query_fasta = self.get_seq(cmd, search_targets=list(query_id))
		out_file = self.blastp(
			wd=self.wd, query_fasta=query_fasta, db=os.path.join(self.db_dir, dmnd)
			)
		out = read_out(out_file)
		hit_id = [hit[0] for hit in out]
		return hit_id, out
	
	def run(self):
		pdb_out = []
		sprot_out = []
		trembl_out = []
		
		# Run pdb
		pdb_cmd = "SELECT query_name, sequence from results_summary"
		pdb_query_fasta = self.get_seq(pdb_cmd)
		all_query_id = [query[0] for query in pdb_query_fasta]
		pdb_out_file = self.blastp(wd=self.wd, query_fasta=pdb_query_fasta, db=os.path.join(self.db_dir, "pdb.dmnd"))
		pdb_out = read_out(pdb_out_file)
		pdb_hit_id = [hit[0] for hit in pdb_out]
		
		# Run swiss prot
		sprot_query_id = set(all_query_id) - set(pdb_hit_id)
		if sprot_query_id:
			sprot_hit_id, sprot_out = self.run_uniprot(query_id=sprot_query_id, dmnd=os.path.join(self.db_dir, "sprot.dmnd"))
			
			# Run trembl
			trembl_query_id = set(sprot_query_id) - set(sprot_hit_id)
			if trembl_query_id:
				trembl_hit_id, trembl_out = self.run_uniprot(query_id=sprot_query_id, dmnd=os.path.join(self.db_dir, "trembl.dmnd"))
		return pdb_out, sprot_out, trembl_out


class UpdateSQL(SqlBuilder, SqlSearch):
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
	wd = r"D:\subject\active\1-qProtein\code\test"
	sql_db = r"D:\subject\active\1-qProtein\code\test\qprotein_results.db"
	db_dir = r"D:\subject\active\1-qProtein\code\db"
	run = RunBlast(wd=wd, sql_db=sql_db, db_dir=db_dir)