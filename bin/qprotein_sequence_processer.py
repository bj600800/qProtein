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


class CazyAnalysis(SqlBuilder):
    def __init__(self, cazy_output, sql_db, table_name, column_definition):
        super(CazyAnalysis, self).__init__(sql_db, table_name, column_definition)

        self.cazy_output = cazy_output
        self.record_items = {'query_name': '', 'ec_number': '', 'cazy_family': ''}

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
            if self.record_counter == 500000:
                self.insert_many(cursor, self.records, 3)
                logger.info(f"Insert 500000 records into {self.table_name}, total records:" + str(self.total_records))

                self.records = []
                self.record_counter = 0
            else:
                line = line.split('\t')
                self.format_records(line)

        if self.record_counter > 0:
            logger.info(f"Insert {self.record_counter} records into {self.table_name}")
            self.insert_many(cursor, self.records, 3)

        logger.info(f"Total records:" + str(self.total_records))

    def run(self):
        logger.info(f"Start to create SQL table: {self.table_name} in SQL file {self.sql_db}")

        logger.info(f"Create SQL table: {self.table_name}")
        cursor = self.create_table(self.table_name, self.sql_db, self.column_definition)

        logger.info(f"Parse file and insert records into {self.table_name}")
        self.parse_cazy(cursor)

        logger.info(f'Creating index for SQL table: {self.table_name}')
        self.create_index_all(cursor)


if __name__ == '__main__':
    cazy_output = r'D:\subject\active\1-qProtein\data\tibet\cazy_overview.txt'
    sql_db = r'C:\Users\bj600\Desktop\qprotein_results.db'
    table_name = 'results_summary'
    column_definition = {'query_name': 'TEXT', 'ec_number': 'TEXT', 'cazy_family': 'TEXT'}
    cazy = CazyAnalysis(cazy_output, sql_db, table_name, column_definition)
    cazy.run()


class MeropsAnalysis(SqlBuilder, SqlSearch):
    def __init__(self, merops_output, sql_db, table_name, column_definition):
        super(MeropsAnalysis, self).__init__(sql_db, table_name, column_definition)
        self.sql_db = sql_db
        self.table_name = table_name

        self.column_definition = column_definition
        self.merops_output = merops_output
        self.record_items = {'query_name': '', 'merops70_family': ''}

        self.records = []
        self.record_counter = 0
        self.total_records = 0

    def parse_items(self, line):
        self.record_items['query_name'] = line[0]
        self.record_items['merops70_family'] = line[1].split('#')[1]

    def format_records(self, line):
        self.parse_items(line)
        self.records.append(tuple(self.record_items.values()))
        self.record_counter += 1
        self.total_records += 1

    def parse_merops(self, cursor, new_column_name, new_column_attr):
        self.add_column(cursor, self.table_name, new_column_name, new_column_attr)
        column = ', '.join(list(self.column_definition.keys()))
        for line in self.read_text_generator(self.merops_output, header=True):
            if self.record_counter == 500000:
                self.insert_many_column(cursor, column, self.records, 2)
                logger.info(f"Insert 500000 records into {self.table_name}, total records:" + str(self.total_records))

                self.records = []
                self.record_counter = 0
            else:
                line = line.split('\t')
                self.format_records(line)

        if self.record_counter > 0:
            logger.info(f"Insert {self.record_counter} records into {self.table_name}")
            self.insert_many_column(cursor, column, self.records, 2)

    def run(self):
        cursor = self.connect_sql(self.sql_db)

        logger.info(f"Parse Merops and insert records into {self.table_name}")
        self.parse_merops(cursor, 'merops70_family', 'TEXT')

        logger.info(f'Creating index for SQL table: {self.table_name}')
        self.create_index_one(cursor, list(self.column_definition.keys())[1])


if __name__ == '__main__':
    merops_output = r'D:\subject\active\1-qProtein\data\tibet\merops_output.tab'
    sql_db = r'C:\Users\bj600\Desktop\qprotein_results.db'
    table_name = 'results_summary'
    column_definition = {'query_name': 'TEXT', 'merops70_family': 'TEXT'}
    merops = MeropsAnalysis(merops_output, sql_db, table_name, column_definition)
    merops.run()
