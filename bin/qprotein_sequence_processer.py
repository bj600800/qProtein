"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: main function to process and manage.
# ------------------------------------------------------------------------------
"""
import sys
sys.path.append("..")

from qprotein.database import annotation_processer
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


def create_sql(summary_sql_path):
    summary_table_name = 'results_summary'
    annotation_processer.CreateSql(summary_table_name, summary_sql_path)
    logger.info(f'Create SQL table {summary_table_name} at {summary_sql_path}')


def insert_cazy(cazy_output, sql_path, summary_table_name):
    column_definition = [('query_name', 'TEXT'), ('ec_number', 'TEXT'), ('cazy_family', 'TEXT')]
    cazy = annotation_processer.CazyAnalysis(cazy_output, sql_path, summary_table_name, column_definition)
    cazy.run()


def insert_merops(merops_output, summary_sql_path, summary_table_name):
    column_definition = [('query_name', 'TEXT'), ('merops_family', 'TEXT')]
    merops = annotation_processer.MeropsAnalysis(merops_output, summary_sql_path, summary_table_name, column_definition)
    merops.run()


def insert_sprot(dmnd_output, summary_sql_path, summary_table_name):
    column_definition = [('query_name', 'TEXT'), ('sprot_acc', 'TEXT'),
                         ('sprot_start', 'TEXT'), ('sprot_end', 'TEXT')]
    dmnd = annotation_processer.SprotDmnd(dmnd_output, summary_sql_path, summary_table_name, column_definition, 'sprot')
    dmnd.run()


def annotate_sprot(uniprot_db, uniprot_table_name, summary_sql_path, summary_table_name):
    column_definition = [('sprot_acc', 'TEXT'), ('sprot_name', 'TEXT'), ('sprot_ec_number', 'TEXT'),
                         ('sprot_go_component', 'TEXT'), ('sprot_go_process', 'TEXT'), ('sprot_go_function', 'TEXT'),
                         ('sprot_interpro', 'TEXT'), ('sprot_pfam', 'TEXT')]

    sprot_dat = annotation_processer.SprotAnnotation(uniprot_db, uniprot_table_name, summary_sql_path,
                                                     summary_table_name, 'sprot_acc', column_definition
                                                     )
    sprot_dat.run()


def insert_trembl(dmnd_output, summary_sql_path, summary_table_name):
    column_definition = [('query_name', 'TEXT'), ('trembl_acc', 'TEXT'),
                         ('trembl_start', 'TEXT'), ('trembl_end', 'TEXT')]

    trembl = annotation_processer.TremblDmnd(dmnd_output, summary_sql_path, summary_table_name, column_definition,
                                             'trembl'
                                             )
    trembl.run()


def annotate_trembl(summary_sql_path, uniprot_db, summary_table_name):
    column_definition = [('trembl_acc', 'TEXT'), ('trembl_name', 'TEXT'), ('trembl_ec_number', 'TEXT'),
                         ('trembl_go_component', 'TEXT'), ('trembl_go_process', 'TEXT'), ('trembl_go_function', 'TEXT'),
                         ('trembl_interpro', 'TEXT'), ('trembl_pfam', 'TEXT')]

    trembl_dat = annotation_processer.TremblAnnotation(summary_sql_path, summary_table_name, uniprot_db,
                                                       column_definition, 'trembl_dat',
                                                       'trembl_acc')
    trembl_dat.run()


def run():
    summary_sql_path = r'D:\subject\active\1-qProtein\data\tibet\qprotein_results.db'
    cazy_output = r'D:\subject\active\1-qProtein\data\tibet\cazy_overview.txt'
    merops_output = r'D:\subject\active\1-qProtein\data\tibet\merops_output.tab'
    sprot_dmnd_output = r'D:\subject\active\1-qProtein\data\tibet\sprot90_90.out'
    summary_table_name = 'results_summary'
    create_sql(summary_sql_path)
    insert_cazy(cazy_output, summary_sql_path, summary_table_name)
    insert_merops(merops_output, summary_sql_path, summary_table_name)
    insert_sprot(sprot_dmnd_output, summary_sql_path, summary_table_name)
    uniprot_db = r'D:\subject\active\1-qProtein\data\tibet\qprotein_db.db'
    uniprot_table_name = 'sprot_dat'
    annotate_sprot(uniprot_db, uniprot_table_name, summary_sql_path, summary_table_name)

    # trembl_dmnd_output = r'D:\subject\active\1-qProtein\data\tibet\trembl_90.tab'
    # insert_trembl(trembl_dmnd_output, summary_sql_path, summary_table_name)
    #
    # annotate_trembl(summary_sql_path, uniprot_db, summary_table_name)


run()
