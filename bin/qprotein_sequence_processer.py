"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: main function to process and manage.
# ------------------------------------------------------------------------------
"""
import os.path
import sys

sys.path.append("..")

from qprotein.database import annotation_processer
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


def create_sql(summary_sql_path):
    annotation_processer.CreateSql(summary_sql_path)
    logger.info(f'Create SQL table results_summary {summary_sql_path}')


def insert_cazy(cazy_output, sql_path, summary_table_name):
    try:
        column_definition = [('query_name', 'TEXT'), ('ec_number', 'TEXT'), ('cazy_family', 'TEXT')]
        cazy = annotation_processer.CazyAnalysis(cazy_output, sql_path, column_definition)
        cazy.run()
    except IOError:
        logger.debug('NO cazy_output file!')


def insert_merops(merops_output, summary_sql_path, summary_table_name):
    try:
        column_definition = [('query_name', 'TEXT'), ('merops_family', 'TEXT')]
        merops = annotation_processer.MeropsAnalysis(merops_output, summary_sql_path, column_definition)
        merops.run()
    except IOError:
        logger.debug('No merops_output file!')


def insert_sprot(dmnd_output, summary_sql_path, summary_table_name):
    try:
        column_definition = [('query_name', 'TEXT'), ('sprot_acc', 'TEXT'),
                             ('sprot_start', 'TEXT'), ('sprot_end', 'TEXT')]
        dmnd = annotation_processer.SprotDmnd(dmnd_output, summary_sql_path, summary_table_name, column_definition, 'sprot')
        dmnd.run()
    except IOError:
        logger.debug('No sprot_dmnd_output file!')


def annotate_sprot(uniprot_db, uniprot_table_name, summary_sql_path, summary_table_name):
    # here create the sprot_acc anain, that is wrong.
    try:
        column_definition = [('sprot_acc', 'TEXT'), ('sprot_name', 'TEXT'), ('sprot_ec_number', 'TEXT'),
                             ('sprot_go_component', 'TEXT'), ('sprot_go_process', 'TEXT'), ('sprot_go_function', 'TEXT'),
                             ('sprot_interpro', 'TEXT'), ('sprot_pfam', 'TEXT')]

        sprot_dat = annotation_processer.SprotAnnotation(uniprot_db, uniprot_table_name, summary_sql_path,
                                                         summary_table_name, 'sprot_acc', column_definition
                                                         )
        sprot_dat.run()
    except Exception as e:
        logger.debug(e)


def insert_trembl(dmnd_output, summary_sql_path, summary_table_name):
    try:
        column_definition = [('query_name', 'TEXT'), ('trembl_acc', 'TEXT'),
                             ('trembl_start', 'TEXT'), ('trembl_end', 'TEXT')]

        trembl = annotation_processer.TremblDmnd(dmnd_output, summary_sql_path, summary_table_name, column_definition,
                                                 'trembl'
                                                 )
        trembl.run()
    except IOError:
        logger.debug('No trembl_dmnd_output file!')


def annotate_trembl(summary_sql_path, uniprot_db, summary_table_name):
    try:
        column_definition = [('trembl_acc', 'TEXT'), ('trembl_name', 'TEXT'), ('trembl_ec_number', 'TEXT'),
                             ('trembl_go_component', 'TEXT'), ('trembl_go_process', 'TEXT'), ('trembl_go_function', 'TEXT'),
                             ('trembl_interpro', 'TEXT'), ('trembl_pfam', 'TEXT')]

        trembl_dat = annotation_processer.TremblAnnotation(summary_sql_path, summary_table_name, uniprot_db,
                                                           column_definition, 'trembl_dat',
                                                           'trembl_acc'
                                                           )
        trembl_dat.run()
    except Exception as e:
        logger.error('Got an exception!', e)


def run():
    root_dir = r'D:\subject\active\1-qProtein\data'
    task_name = 'tibet'
    work_dir = os.path.join(root_dir, task_name)
    if not os.path.exists(work_dir):
        logger.info('Create dir', work_dir)
        os.mkdir(work_dir)
    summary_sql_path = os.path.join(work_dir, 'qprotein_results.db')
    cazy_output = os.path.join(work_dir, 'cazy_overview.txt')
    merops_output = os.path.join(work_dir, 'merops_output.tab')
    sprot_dmnd_output = os.path.join(work_dir, 'sprot90_90.out')
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