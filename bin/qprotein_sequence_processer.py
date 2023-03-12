"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: perform sequence annotation and sql builder for each database
# ------------------------------------------------------------------------------
"""

from qprotein.database.sqlite3_builder import SqlBuilder
from qprotein.database.sqlite3_searcher import SqlSearch
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


class CreateSql(SqlBuilder):
    def __init__(self, table_name, sql_db):
        super(SqlBuilder, self).__init__()
        self.table_name = table_name
        self.sql_db = sql_db

        self.create()

    def create(self):
        sql = f"CREATE TABLE {self.table_name} (query_name TEXT PRIMARY KEY)"
        self.create_table(self.table_name, self.sql_db, sql)


table_name = 'results_summary'
sql_db = r'C:\Users\bj600\Desktop\qprotein_results.db'

results_sql = CreateSql(table_name, sql_db)
logger.info(f'Create SQL table {table_name} at {sql_db}')


class CazyAnalysis(SqlBuilder):
    def __init__(self, cazy_output, sql_db, table_name, column_definition):
        super(CazyAnalysis, self).__init__(sql_db, table_name, column_definition)

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
        for line in self.read_text_generator(self.cazy_output, header=True):
            line = line.split('\t')

            if self.record_counter == 500000:
                self.insert_many(cursor, self.records, 3)
                logger.info(f"Insert 500000 records into {self.table_name}, total records:" + str(self.total_records))

                self.records = []
                self.record_counter = 0

            else:
                self.format_records(line)

        if self.record_counter > 0:
            logger.info(f"Insert {self.record_counter} records into {self.table_name}")
            self.insert_many(cursor, self.records, 3)

        logger.info(f"Total records:" + str(self.total_records))

    def run(self):
        columns = [i[0] for i in self.column_definition]
        logger.info(f"Connect SQL table: {self.table_name} in {self.sql_db}")
        cursor = self.connect_sql(self.sql_db)

        logger.info(f"Add columns {', '.join(columns)} into table {self.table_name}")
        self.add_column(cursor, self.table_name, self.column_definition)

        logger.info(f"!! Parse cazy output and insert records")
        self.parse_cazy(cursor)

        logger.info(f'Create index cazy_idx for SQL table: {self.table_name}')
        self.create_index(cursor, columns)


if __name__ == '__main__':
    cazy_output = r'D:\subject\active\1-qProtein\data\tibet\cazy_overview.txt'

    table_name = 'results_summary'
    column_definition = [('query_name', 'TEXT'), ('ec_number', 'TEXT'), ('cazy_family', 'TEXT')]
    cazy = CazyAnalysis(cazy_output, sql_db, table_name, column_definition)
    cazy.run()


class MeropsAnalysis(SqlBuilder, SqlSearch):
    def __init__(self, merops_output, sql_db, table_name, column_definition):
        super(MeropsAnalysis, self).__init__(sql_db, table_name, column_definition)
        self.sql_db = sql_db
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

        logger.info(f"Connect SQL table: {self.table_name} in {self.sql_db}")
        cursor = self.connect_sql(self.sql_db)

        logger.info(f"Add columns {', '.join(columns)} into table {self.table_name}")
        self.add_column(cursor, self.table_name, self.column_definition)

        logger.info(f"!! Parse Merops output and insert records")
        self.parse_merops(cursor, columns)

        logger.info(f'Create index merops_idx for SQL table: {self.table_name}')
        self.create_index(cursor, columns)


if __name__ == '__main__':
    merops_output = r'D:\subject\active\1-qProtein\data\tibet\merops_output.tab'
    sql_db = r'C:\Users\bj600\Desktop\qprotein_results.db'
    table_name = 'results_summary'
    column_definition = [('query_name', 'TEXT'), ('merops_family', 'TEXT')]
    merops = MeropsAnalysis(merops_output, sql_db, table_name, column_definition)
    merops.run()


class SprotDmnd(SqlBuilder, SqlSearch):
    def __init__(self, dmnd_output, sql_db, table_name, column_definition, idx_name_prefix):
        super().__init__(sql_db, table_name, column_definition)

        self.dmnd_output = dmnd_output
        self.sql_db = sql_db
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
        if self.record_items['sprot_acc'] == 'P28074':
            print(line)
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

        logger.info(f"Connect SQL table: {self.table_name} in {self.sql_db}")
        cursor = self.connect_sql(self.sql_db)

        logger.info(f"Add columns {', '.join(columns)} into table {self.table_name}")
        self.add_column(cursor, self.table_name, self.column_definition)

        logger.info(f"!! Parse {self.idx_name_prefix} and insert records into {self.table_name}")
        self.parse_dmnd(cursor, columns)

        logger.info(f'Create index for SQL table: {self.table_name}')
        self.create_index(cursor, columns)


if __name__ == '__main__':
    dmnd_output = r'D:\subject\active\1-qProtein\data\tibet\sprot_output_70.tab'
    sql_db = r'C:\Users\bj600\Desktop\qprotein_results.db'
    table_name = 'results_summary'
    column_definition = [('query_name', 'TEXT'), ('sprot_acc', 'TEXT'),
                         ('sprot_start', 'TEXT'), ('sprot_end', 'TEXT')]

    dmnd = SprotDmnd(dmnd_output, sql_db, table_name, column_definition, 'sprot')
    dmnd.run()


class SprotAnnotation(SqlBuilder, SqlSearch):
    def __init__(self, sql_db, table_name, uniprot_db, column_definition, uniprot_table, target_column):
        super().__init__(sql_db, table_name, column_definition)

        self.sql_db = sql_db
        self.uniprot_db = uniprot_db
        self.table_name = table_name
        self.column_definition = column_definition
        self.uniprot_table = uniprot_table
        self.target_column = target_column

        self.record_items = {column_name[0]: '' for column_name in self.column_definition}
        self.record_counter = 0
        self.total_records = 0

    def get_candidates(self):
        cursor = self.connect_sql(self.sql_db)
        sql_cmd = f"SELECT {self.target_column} FROM results_summary WHERE {self.target_column} != ''"
        acc = self.fetch_results(cursor, sql_cmd)
        return acc

    def get_uniprot_dat(self, uniprot_table_name):
        acc = self.get_candidates()
        cursor = self.connect_sql(self.uniprot_db)
        # prepare
        batch_size = 100000
        uniprot_dat = []
        total_query = len(acc)
        for i in range(0, total_query, batch_size):
            batch_acc = acc[i:i + batch_size]
            sql_cmd = "SELECT * FROM {} WHERE accession IN ({})".format(uniprot_table_name,
                                                                        ', '.join(['?'] * len(batch_acc))
                                                                        )
            cursor.execute(sql_cmd, batch_acc)
            return_dat = cursor.fetchall()
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
        logger.info(f"Connect SQL table: {self.table_name} in {self.sql_db}")
        columns = [i[0] for i in self.column_definition]

        cursor = self.connect_sql(self.sql_db)

        logger.info(f"Add columns {', '.join(columns)} into table {self.table_name}")
        self.add_column(cursor, self.table_name, self.column_definition)

        logger.info(f"!! Parse sprot: update and insert records into {self.table_name}")
        self.parse_uniprot(cursor, columns)

        logger.info(f'Create index for SQL table: {self.table_name}')
        self.create_index(cursor, columns)


if __name__ == '__main__':
    sql_db = r'C:\Users\bj600\Desktop\qprotein_results.db'
    uniprot_db = r'C:\Users\bj600\Desktop\qprotein_db.db'
    table_name = 'results_summary'
    uniprot_table = 'sprot_dat'
    target_column = 'sprot_acc'
    column_definition = [('sprot_acc', 'TEXT'), ('sprot_name', 'TEXT'), ('sprot_ec_number', 'TEXT'),
                         ('sprot_go_component', 'TEXT'), ('sprot_go_process', 'TEXT'), ('sprot_go_function', 'TEXT'),
                         ('sprot_interpro', 'TEXT'), ('sprot_pfam', 'TEXT')]

    sprot_dat = SprotAnnotation(sql_db, table_name, uniprot_db, column_definition, uniprot_table, target_column)
    sprot_dat.run()


class TremblDmnd(SprotDmnd):
    def __init__(self, dmnd_output, sql_db, table_name, column_definition, idx_name_prefix):
        super().__init__(dmnd_output, sql_db, table_name, column_definition, idx_name_prefix)
        self.dmnd_output = dmnd_output
        self.sql_db = sql_db
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


if __name__ == '__main__':
    dmnd_output = r'D:\subject\active\1-qProtein\data\tibet\trembl_output_70.tab'
    sql_db = r'C:\Users\bj600\Desktop\qprotein_results.db'
    table_name = 'results_summary'
    column_definition = [('query_name', 'TEXT'), ('trembl_acc', 'TEXT'),
                         ('trembl_start', 'TEXT'), ('trembl_end', 'TEXT')]

    trembl = TremblDmnd(dmnd_output, sql_db, table_name, column_definition, 'trembl')
    trembl.run()


class TremblAnnotation(SprotAnnotation):
    def __init__(self, sql_db, table_name, uniprot_db, column_definition, uniprot_table, target_column):
        super().__init__(sql_db, table_name, uniprot_db, column_definition, uniprot_table, target_column)

        self.sql_db = sql_db
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
        logger.info(f"Connect SQL table: {self.table_name} in {self.sql_db}")
        cursor = self.connect_sql(self.sql_db)

        logger.info(f"Add columns {', '.join(columns)} into table {self.table_name}")
        self.add_column(cursor, self.table_name, self.column_definition)

        logger.info(f"!! Parse Trembl update and insert records into {self.table_name}")
        self.parse_uniprot(cursor, columns)

        logger.info(f'Create index for SQL table: {self.table_name}')

        self.create_index(cursor, columns)


if __name__ == '__main__':
    sql_db = r'C:\Users\bj600\Desktop\qprotein_results.db'
    uniprot_db = r'C:\Users\bj600\Desktop\qprotein_db.db'
    table_name = 'results_summary'
    uniprot_table = 'trembl_dat'
    target_column = 'trembl_acc'
    column_definition = [('trembl_acc', 'TEXT'), ('trembl_name', 'TEXT'), ('trembl_ec_number', 'TEXT'),
                         ('trembl_go_component', 'TEXT'), ('trembl_go_process', 'TEXT'), ('trembl_go_function', 'TEXT'),
                         ('trembl_interpro', 'TEXT'), ('trembl_pfam', 'TEXT')]

    trembl_dat = TremblAnnotation(sql_db, table_name, uniprot_db, column_definition, uniprot_table, target_column)
    trembl_dat.run()
