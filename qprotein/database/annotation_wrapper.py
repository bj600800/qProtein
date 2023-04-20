"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: process sequence annotation outputs and build sql database
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
    def __init__(self, dmnd_output, summary_sql_path):
        super().__init__(sql_db=summary_sql_path, table_name="results_summary",
                         column_definition=[('query_name', 'TEXT'), ('sprot_acc', 'TEXT'),
                                            ('sprot_ident', 'FLOAT'), ('sprot_cover', 'FLOAT'),
                                            ('sprot_match_length', 'INT'), ('sprot_query_start', 'INT'),
                                            ('sprot_query_end', 'INT'), ('sprot_subject_start', 'INT'),
                                            ('sprot_subject_end', 'INT'), ('sprot_evalue', 'FLOAT')])
        self.dmnd_output = dmnd_output
        self.summary_sql_path = summary_sql_path
        self.idx_name_prefix = 'sprot'

        self.record_items = {column_name[0]: '' for column_name in self.column_definition}
        self.records = []
        self.record_counter = 0
        self.total_records = 0

    def parse_items_from_dmnd(self, line):
        self.record_items['query_name'] = line[0]
        self.record_items['sprot_acc'] = line[1].split('|')[1]
        self.record_items['sprot_ident'] = line[2]
        self.record_items['sprot_match_length'] = line[3]
        self.record_items['sprot_query_start'] = line[4]
        self.record_items['sprot_query_end'] = line[5]
        self.record_items['sprot_subject_start'] = line[6]
        self.record_items['sprot_subject_end'] = line[7]
        self.record_items['sprot_cover'] = line[8]
        self.record_items['sprot_evalue'] = line[9].rstrip()

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
                logger.info(f"Insert 500000 records. Total records:" + str(self.total_records))
                self.insert_update_columns(cursor=cursor, table_name="results_summary", columns=columns,
                                           records=self.records
                                           )
                self.records = []
                self.record_counter = 0
            else:
                self.format_dmnd(line)

        if self.record_counter > 0:
            logger.info(f"Insert {self.record_counter} records")
            self.insert_update_columns(cursor=cursor, table_name="results_summary", columns=columns,
                                       records=self.records
                                       )

        logger.info(f"Total records:" + str(self.total_records))

    def run(self):
        logger.info(f"Parse {self.idx_name_prefix} diamond outputs and insert records")
        columns = [i[0] for i in self.column_definition]
        cursor = self.connect_sql(sql_db=self.summary_sql_path)
        self.add_column(cursor=cursor, table_name="results_summary", column_definition=self.column_definition)
        self.parse_dmnd(cursor=cursor, columns=columns)
        self.create_index(cursor=cursor, columns=columns)


class SprotAnnotation(SqlBuilder, SqlSearch):
    def __init__(self, uniprot_db, summary_sql_path):
        super().__init__(sql_db=summary_sql_path, table_name="results_summary",
                         column_definition=[('sprot_acc', 'TEXT'), ('sprot_name', 'TEXT'), ('sprot_ec_number', 'TEXT'),
                                            ('sprot_go_component', 'TEXT'), ('sprot_go_process', 'TEXT'),
                                            ('sprot_go_function', 'TEXT'), ('sprot_interpro', 'TEXT'),
                                            ('sprot_pfam', 'TEXT')])

        self.summary_sql_path = summary_sql_path
        self.uniprot_db = uniprot_db
        self.uniprot_table_name = "sprot_dat"
        self.target_column = "sprot_acc"

        self.record_items = {column_name[0]: '' for column_name in self.column_definition}
        self.record_counter = 0
        self.total_records = 0

    def get_candidates(self):
        cursor = self.connect_sql(self.summary_sql_path)
        sql_cmd = f"SELECT {self.target_column} FROM results_summary WHERE {self.target_column} != ''"
        acc = [i[0] for i in self.fetch_results(cursor=cursor, sql_cmd=sql_cmd)]
        return acc

    def get_uniprot_dat(self, uniprot_table_name):
        acc = self.get_candidates()
        uniprot_cursor = self.connect_sql(sql_db=self.uniprot_db)
        # prepare
        batch_size = 100000
        uniprot_dat = []
        total_query = len(acc)
        for i in range(0, total_query, batch_size):
            batch_acc = acc[i:i + batch_size]
            sql_cmd = "SELECT * FROM {} WHERE accession IN ({})" \
                .format(uniprot_table_name, ', '.join(['?'] * len(batch_acc)))
            return_dat = SqlSearch.fetch_many_results(cursor=uniprot_cursor, sql_cmd=sql_cmd, search_targets=batch_acc)
            uniprot_dat.extend(return_dat)
        return uniprot_dat

    def parse_items_from_db(self, dat):
        self.record_items['sprot_acc'] = dat[0]
        self.record_items['sprot_name'] = dat[1]
        self.record_items['sprot_ec_number'] = dat[2]
        self.record_items['sprot_go_component'] = dat[3]
        self.record_items['sprot_go_process'] = dat[4]
        self.record_items['sprot_go_function'] = dat[5]
        self.record_items['sprot_interpro'] = dat[6]
        self.record_items['sprot_pfam'] = dat[7]

    def format_records(self, dat):
        self.parse_items_from_db(dat)
        self.records.append(tuple(self.record_items.values()))
        self.record_counter += 1
        self.total_records += 1

    def parse_uniprot(self, cursor, columns):
        uniprot_dat = self.get_uniprot_dat(uniprot_table_name=self.uniprot_table_name)
        for i in uniprot_dat:
            if self.record_counter == 500000:
                self.update_many_columns(cursor=cursor, table_name="results_summary",
                                         columns=columns, records=self.records
                                         )
                logger.info(f"Insert 500000 records. Total records:" + str(self.total_records))

                self.records = []
                self.record_counter = 0
            else:
                self.format_records(i)
        if self.record_counter > 0:
            self.update_many_columns(cursor=cursor, table_name="results_summary", columns=columns, records=self.records)
            logger.info(f"Update {self.record_counter} records")
        logger.info(f"Total records:" + str(self.total_records))

    def run(self):
        logger.info(f"Annotate based on sprot accession")
        columns = [i[0] for i in self.column_definition]
        cursor = self.connect_sql(sql_db=self.summary_sql_path)
        # annotation should not add accession column
        self.add_column(cursor=cursor, table_name="results_summary", column_definition=self.column_definition[1:])
        self.parse_uniprot(cursor=cursor, columns=columns)
        self.create_index(cursor=cursor, columns=columns)


#
class TremblDmnd(SprotDmnd):
    def __init__(self, dmnd_output, summary_sql_path):
        super().__init__(dmnd_output=dmnd_output, summary_sql_path=summary_sql_path,

                         )
        self.dmnd_output = dmnd_output
        self.summary_sql_path = summary_sql_path
        self.column_definition = [('query_name', 'TEXT'), ('trembl_acc', 'TEXT'),
                                  ('trembl_ident', 'FLOAT'), ('trembl_cover', 'FLOAT'),
                                  ('trembl_match_length', 'INT'), ('trembl_query_start', 'INT'),
                                  ('trembl_query_end', 'INT'), ('trembl_subject_start', 'INT'),
                                  ('trembl_subject_end', 'INT'), ('trembl_evalue', 'FLOAT')]

        self.idx_name_prefix = 'trembl'
        self.record_items = {column_name[0]: '' for column_name in self.column_definition}
        self.records = []
        self.record_counter = 0
        self.total_records = 0

    def parse_items_from_dmnd(self, line):
        self.record_items['query_name'] = line[0]
        self.record_items['trembl_acc'] = line[1].split('|')[1]
        self.record_items['trembl_ident'] = line[2]
        self.record_items['trembl_match_length'] = line[3]
        self.record_items['trembl_query_start'] = line[4]
        self.record_items['trembl_query_end'] = line[5]
        self.record_items['trembl_subject_start'] = line[6]
        self.record_items['trembl_subject_end'] = line[7]
        self.record_items['trembl_cover'] = line[8]
        self.record_items['trembl_evalue'] = line[9]


class TremblAnnotation(SprotAnnotation):
    def __init__(self, summary_sql_path, uniprot_db):
        super().__init__(summary_sql_path=summary_sql_path, uniprot_db=uniprot_db)

        self.summary_sql_path = summary_sql_path
        self.uniprot_db = uniprot_db
        self.column_definition = [('trembl_acc', 'TEXT'), ('trembl_name', 'TEXT'), ('trembl_ec_number', 'TEXT'),
                                  ('trembl_go_component', 'TEXT'), ('trembl_go_process', 'TEXT'),
                                  ('trembl_go_function', 'TEXT'), ('trembl_interpro', 'TEXT'), ('trembl_pfam', 'TEXT')]
        self.uniprot_table_name = "trembl_dat"
        self.target_column = "trembl_acc"

        self.record_items = {column_name[0]: '' for column_name in self.column_definition}
        self.record_counter = 0
        self.total_records = 0

    def parse_items_from_db(self, dat):
        self.record_items['trembl_acc'] = dat[0]
        self.record_items['trembl_name'] = dat[1]
        self.record_items['trembl_ec_number'] = dat[2]
        self.record_items['trembl_go_component'] = dat[3]
        self.record_items['trembl_go_process'] = dat[4]
        self.record_items['trembl_go_function'] = dat[5]
        self.record_items['trembl_interpro'] = dat[6]
        self.record_items['trembl_pfam'] = dat[7]

    def run(self):
        logger.info(f"Annotate based on trembl accession")
        columns = [i[0] for i in self.column_definition]
        cursor = self.connect_sql(sql_db=self.summary_sql_path)
        self.add_column(cursor=cursor, table_name="results_summary", column_definition=self.column_definition[1:])
        self.parse_uniprot(cursor=cursor, columns=columns)
        self.create_index(cursor=cursor, columns=columns)


class CazyAnalysis(SqlBuilder):
    def __init__(self, cazy_output, sql_path):
        super(CazyAnalysis, self).__init__(sql_db=sql_path, table_name='results_summary',
                                           column_definition=[('query_name', 'TEXT'), ('ec_number', 'TEXT'),
                                                              ('cazy_family', 'TEXT')])
        self.cazy_output = cazy_output
        self.record_items = {column_name[0]: '' for column_name in self.column_definition}

        self.records = []
        self.record_counter = 0
        self.total_records = 0

    def parse_items(self, line):
        self.record_items['query_name'] = line[0]
        if line[1] == '-':
            self.record_items['ec_number'] = ''
        else:
            self.record_items['ec_number'] = line[1]
        for i in line[2:5]:
            if i != '-':
                self.record_items['cazy_family'] = i

    def format_records(self, line):
        self.parse_items(line)
        self.records.append(tuple(self.record_items.values()))
        self.record_counter += 1
        self.total_records += 1

    def parse_cazy(self, cursor, columns):
        for line in self.read_text_generator(filename=self.cazy_output, header=True):
            line = line.split('\t')
            if self.record_counter == 500000:
                self.insert_update_columns(cursor=cursor, table_name="results_summary", columns=columns, records=self.records)
                logger.info(f"Insert 500000 records. Total records:" + str(self.total_records))

                self.records = []
                self.record_counter = 0
            else:
                self.format_records(line)

        if self.record_counter > 0:
            logger.info(f"Insert {self.record_counter} records")
            self.insert_update_columns(cursor=cursor, table_name="results_summary", columns=columns,
                                       records=self.records
                                       )

        logger.info(f"Total records:" + str(self.total_records))

    def run(self):
        logger.info(f"Parse Cazy output and insert records")
        columns = [i[0] for i in self.column_definition]
        cursor = self.connect_sql(sql_db=self.sql_path)
        self.add_column(cursor=cursor, table_name='results_summary', column_definition=self.column_definition)
        self.parse_cazy(cursor=cursor, columns=columns)
        self.create_index(cursor=cursor, columns=columns)


class MeropsAnalysis(SqlBuilder, SqlSearch):
    def __init__(self, merops_output, summary_sql_path):
        super(MeropsAnalysis, self).__init__(sql_db=summary_sql_path, table_name='results_summary',
                                             column_definition=[('query_name', 'TEXT'), ('merops_family', 'TEXT')])
        self.summary_sql_path = summary_sql_path
        self.merops_output = merops_output
        self.record_items = {column_name[0]: '' for column_name in self.column_definition}

        self.records = []
        self.record_counter = 0
        self.total_records = 0

    def parse_items(self, line):
        self.record_items['query_name'] = line[0]
        self.record_items['merops_family'] = line[1].split('#')[1]

    def format_records(self, line):
        self.parse_items(line)
        self.records.append(tuple(self.record_items.values()))
        self.record_counter += 1
        self.total_records += 1

    def parse_merops(self, cursor, columns):
        for line in self.read_text_generator(filename=self.merops_output, header=True):
            line = line.split('\t')

            if self.record_counter == 500000:
                self.insert_update_columns(cursor=cursor, table_name="results_summary",
                                           columns=columns, records=self.records)
                logger.info(f"Insert 500000 records. Total records:" + str(self.total_records))

                self.records = []
                self.record_counter = 0
            else:
                self.format_records(line)

        if self.record_counter > 0:
            logger.info(f"Insert {self.record_counter} records")
            self.insert_update_columns(cursor=cursor, table_name="results_summary", columns=columns,
                                       records=self.records)

        logger.info(f"Total records:" + str(self.total_records))

    def run(self):
        logger.info(f"Parse Merops output and insert records")
        columns = [i[0] for i in self.column_definition]
        cursor = self.connect_sql(sql_db=self.summary_sql_path)
        self.add_column(cursor=cursor, table_name="results_summary", column_definition=self.column_definition)
        self.parse_merops(cursor=cursor, columns=columns)
        self.create_index(cursor=cursor, columns=columns)


class QueryAnalysis(SqlBuilder, SqlSearch):
    def __init__(self, task_name, summary_sql_path, fasta_sql_path):
        super().__init__(sql_db=summary_sql_path, table_name="results_summary",
                         column_definition=[('query_name', 'TEXT'), ('seq_length', 'INT')])
        self.summary_sql_path = summary_sql_path
        self.fasta_sql_path = fasta_sql_path
        self.record_items = {column_name[0]: '' for column_name in self.column_definition}
        self.task_name = task_name

    def get_query_name(self):
        summary_cursor = self.connect_sql(self.summary_sql_path)
        sql_cmd = "SELECT query_name FROM results_summary WHERE sprot_acc != '' or trembl_acc != ''"
        query_name = [i[0] for i in SqlSearch.fetch_results(summary_cursor, sql_cmd)]
        return query_name

    def get_query_length(self):
        fasta_cursor = self.connect_sql(self.fasta_sql_path)
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
                self.update_many_columns(cursor=cursor, table_name="results_summary",
                                         columns=columns, records=self.records)
                logger.info(f"Insert 500000 records. Total records:" + str(self.total_records))

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
        # annotation should not add accession column
        self.add_column(cursor=cursor, table_name="results_summary",
                        column_definition=self.column_definition[1:], specific_column='query_name')
        self.parse_uniprot(cursor=cursor, columns=columns)
        self.create_index(cursor=cursor, columns=columns)
