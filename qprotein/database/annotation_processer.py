"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
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


class CazyAnalysis(SqlBuilder):
    def __init__(self, cazy_output, sql_path, column_definition):
        super(CazyAnalysis, self).__init__(sql_db=sql_path, table_name='results_summary',
                                           column_definition=column_definition)
        self.cazy_output = cazy_output
        self.column_definition = column_definition
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

    def parse_cazy(self, cursor):
        for line in self.read_text_generator(filename=self.cazy_output, header=True):
            line = line.split('\t')

            if self.record_counter == 500000:
                self.insert_many(cursor=cursor, records=self.records, values_num=3)
                logger.info(f"Insert 500000 records. Total records:" + str(self.total_records))

                self.records = []
                self.record_counter = 0

            else:
                self.format_records(line)

        if self.record_counter > 0:
            logger.info(f"Insert {self.record_counter} records")
            self.insert_many(cursor=cursor, records=self.records, values_num=3)

        logger.info(f"Total records:" + str(self.total_records))

    def run(self):
        columns = [i[0] for i in self.column_definition]
        cursor = self.connect_sql(sql_db=self.sql_path)

        logger.info(f"Add columns {', '.join(columns)} into the table results_summary")
        self.add_column(cursor=cursor, table_name='results_summary', column_definition=self.column_definition)

        logger.info(f"Parse cazy output and insert records")
        self.parse_cazy(cursor=cursor)

        logger.info(f'Create index for cazy')
        self.create_index(cursor=cursor, columns=columns)


class MeropsAnalysis(SqlBuilder, SqlSearch):
    def __init__(self, merops_output, summary_sql_path, table_name, column_definition):
        super(MeropsAnalysis, self).__init__(summary_sql_path, table_name, column_definition)
        self.summary_sql_path = summary_sql_path
        self.table_name = table_name

        self.column_definition = column_definition
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
        for line in self.read_text_generator(self.merops_output, header=True):
            line = line.split('\t')

            if self.record_counter == 500000:
                self.insert_update_columns(cursor, self.table_name, columns, self.records)
                logger.info(f"Insert 500000 records into {self.table_name}, total records:" + str(self.total_records))

                self.records = []
                self.record_counter = 0
            else:
                self.format_records(line)

        if self.record_counter > 0:
            logger.info(f"Insert {self.record_counter} records into {self.table_name}")
            self.insert_update_columns(cursor, self.table_name, columns, self.records)

        logger.info(f"Total records:" + str(self.total_records))

    def run(self):
        columns = [i[0] for i in self.column_definition]

        logger.info(f"Connect SQL table: {self.table_name} in {self.summary_sql_path}")
        cursor = self.connect_sql(self.summary_sql_path)

        logger.info(f"Add columns {', '.join(columns)} into table {self.table_name}")
        self.add_column(cursor, self.table_name, self.column_definition)

        logger.info(f"!! Parse Merops output and insert records")
        self.parse_merops(cursor, columns)

        logger.info(f'Create index merops_idx for SQL table: {self.table_name}')
        self.create_index(cursor, columns)


class SprotDmnd(SqlBuilder, SqlSearch):
    def __init__(self, dmnd_output, summary_sql_path, table_name, column_definition, idx_name_prefix):
        super().__init__(summary_sql_path, table_name, column_definition)

        self.dmnd_output = dmnd_output
        self.summary_sql_path = summary_sql_path
        self.table_name = table_name
        self.column_definition = column_definition
        self.idx_name_prefix = idx_name_prefix

        self.record_items = {column_name[0]: '' for column_name in self.column_definition}
        self.records = []
        self.record_counter = 0
        self.total_records = 0

    def parse_items_from_dmnd(self, line):
        self.record_items['query_name'] = line[0]
        self.record_items['sprot_acc'] = line[1].split('|')[1]
        self.record_items['sprot_start'] = line[8]
        self.record_items['sprot_end'] = line[9]

    def format_dmnd(self, line):
        # prepare
        self.parse_items_from_dmnd(line)
        self.records.append(tuple(self.record_items.values()))
        self.record_counter += 1
        self.total_records += 1

    def parse_dmnd(self, cursor, columns):
        for line in self.read_text_generator(self.dmnd_output, header=False):
            line = line.split('\t')
            if self.record_counter == 500000:
                logger.info(f"Insert 500000 records into {self.table_name}, total records:" + str(self.total_records))
                self.insert_update_columns(cursor, self.table_name, columns, self.records)

                self.records = []
                self.record_counter = 0
            else:
                self.format_dmnd(line)

        if self.record_counter > 0:
            logger.info(f"Insert {self.record_counter} records into {self.table_name}")
            self.insert_update_columns(cursor, self.table_name, columns, self.records)

        logger.info(f"Total records:" + str(self.total_records))

    def run(self):
        columns = [i[0] for i in self.column_definition]

        logger.info(f"Connect SQL table: {self.table_name} in {self.summary_sql_path}")
        cursor = self.connect_sql(self.summary_sql_path)

        logger.info(f"Add columns {', '.join(columns)} into table {self.table_name}")
        self.add_column(cursor, self.table_name, self.column_definition)

        logger.info(f"!! Parse {self.idx_name_prefix} and insert records into {self.table_name}")
        self.parse_dmnd(cursor, columns)

        logger.info(f'Create index for SQL table: {self.table_name}')
        self.create_index(cursor, columns)


class SprotAnnotation(SqlBuilder, SqlSearch):
    def __init__(self, uniprot_db, uniprot_table_name, summary_sql_path, summary_table_name, search_column,
                 column_definition
                 ):
        super().__init__(sql_db=summary_sql_path, table_name=summary_table_name, column_definition=column_definition)

        self.summary_sql_path = summary_sql_path
        self.uniprot_db = uniprot_db
        self.table_name = summary_table_name
        self.column_definition = column_definition
        self.uniprot_table = uniprot_table_name
        self.target_column = search_column

        self.record_items = {column_name[0]: '' for column_name in self.column_definition}
        self.record_counter = 0
        self.total_records = 0

    def get_candidates(self):
        cursor = self.connect_sql(self.summary_sql_path)
        sql_cmd = f"SELECT {self.target_column} FROM results_summary WHERE {self.target_column} != ''"
        acc = [i[0] for i in self.fetch_results(cursor, sql_cmd)]
        return acc

    def get_uniprot_dat(self, uniprot_table_name):
        acc = self.get_candidates()
        uniprot_cursor = self.connect_sql(self.uniprot_db)
        # prepare
        batch_size = 100000
        uniprot_dat = []
        total_query = len(acc)
        for i in range(0, total_query, batch_size):
            batch_acc = acc[i:i + batch_size]
            sql_cmd = "SELECT * FROM {} WHERE accession IN ({})" \
                .format(uniprot_table_name, ', '.join(['?'] * len(batch_acc)))
            uniprot_cursor.execute(sql_cmd, batch_acc)
            return_dat = uniprot_cursor.fetchall()
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
        uniprot_dat = self.get_uniprot_dat(self.uniprot_table)
        for i in uniprot_dat:
            if self.record_counter == 500000:
                self.update_many_columns(cursor, self.table_name, columns, self.records)
                logger.info(f"Insert 500000 records into {self.table_name}, total records:" + str(self.total_records))

                self.records = []
                self.record_counter = 0
            else:
                self.format_records(i)
        if self.record_counter > 0:
            logger.info(f"Update {self.record_counter} records into {self.table_name}")
            self.update_many_columns(cursor, self.table_name, columns, self.records)

    def run(self):
        logger.info(f"Connect SQL table: {self.table_name} in {self.summary_sql_path}")
        columns = [i[0] for i in self.column_definition]

        cursor = self.connect_sql(self.summary_sql_path)

        logger.info(f"Add columns {', '.join(columns)} into table {self.table_name}")
        # annotation should not add accession column
        self.add_column(cursor, self.table_name, self.column_definition[1:])

        logger.info(f"!! Parse sprot: update and insert records into {self.table_name}")
        self.parse_uniprot(cursor, columns)

        logger.info(f'Create index for SQL table: {self.table_name}')
        self.create_index(cursor, columns)


#
class TremblDmnd(SprotDmnd):
    def __init__(self, dmnd_output, summary_sql_path, table_name, column_definition, idx_name_prefix):
        super().__init__(dmnd_output, summary_sql_path, table_name, column_definition, idx_name_prefix)
        self.dmnd_output = dmnd_output
        self.summary_sql_path = summary_sql_path
        self.table_name = table_name
        self.column_definition = column_definition
        self.idx_name_prefix = idx_name_prefix

        self.record_items = {column_name[0]: '' for column_name in self.column_definition}
        self.records = []
        self.record_counter = 0
        self.total_records = 0

    def parse_items_from_dmnd(self, line):
        self.record_items['query_name'] = line[0]
        self.record_items['trembl_acc'] = line[1].split('|')[1]
        self.record_items['trembl_start'] = line[8]
        self.record_items['trembl_end'] = line[9]


class TremblAnnotation(SprotAnnotation):
    def __init__(self, summary_sql_path, table_name, uniprot_db, column_definition, uniprot_table, target_column):
        super().__init__(summary_sql_path, table_name, uniprot_db, column_definition, uniprot_table, target_column)

        self.summary_sql_path = summary_sql_path
        self.uniprot_db = uniprot_db
        self.table_name = table_name
        self.column_definition = column_definition
        self.uniprot_table = uniprot_table
        self.target_column = target_column

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
        columns = [i[0] for i in self.column_definition]
        logger.info(f"Connect SQL table: {self.table_name} in {self.summary_sql_path}")
        cursor = self.connect_sql(self.summary_sql_path)

        logger.info(f"Add columns {', '.join(columns)} into table {self.table_name}")
        self.add_column(cursor, self.table_name, self.column_definition[1:])

        logger.info(f"!! Parse Trembl update and insert records into {self.table_name}")
        self.parse_uniprot(cursor, columns)

        logger.info(f'Create index for SQL table: {self.table_name}')

        self.create_index(cursor, columns)
