"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: main function to process and manage.
# ------------------------------------------------------------------------------
"""
import os.path
import sys
import argparse
sys.path.append("..")

from qprotein.database import annotation_wrapper
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)

parser = argparse.ArgumentParser(description='Crawl query structure with multiprocessing.')
parser.add_argument('--work_dir', required=True, help='All task dirs should be here')
parser.add_argument('--task_name', required=True, help='Name the task')
args = parser.parse_args()


def create_sql(summary_sql_path):
    annotation_wrapper.CreateSql(summary_sql_path)
    logger.info(f"Create SQL table results_summary {summary_sql_path}")


def insert_sprot(dmnd_output, summary_sql_path):
    try:
        dmnd = annotation_wrapper.SprotDmnd(dmnd_output, summary_sql_path)
        dmnd.run()
    except IOError:
        logger.debug("No sprot_dmnd_output file!")


def annotate_sprot(uniprot_db, summary_sql_path):
    try:
        sprot_dat = annotation_wrapper.SprotAnnotation(uniprot_db=uniprot_db, summary_sql_path=summary_sql_path)
        sprot_dat.run()
    except Exception as e:
        logger.debug(e)


def insert_trembl(dmnd_output, summary_sql_path):
    try:
        trembl = annotation_wrapper.TremblDmnd(dmnd_output, summary_sql_path)
        trembl.run()
    except IOError:
        logger.debug("No trembl_dmnd_output file!")


def annotate_trembl(summary_sql_path, uniprot_db):
    try:
        trembl_dat = annotation_wrapper.TremblAnnotation(summary_sql_path=summary_sql_path, uniprot_db=uniprot_db)
        trembl_dat.run()
    except Exception as e:
        logger.error("Got an exception!", e)


def insert_cazy(cazy_output, sql_path):
    try:
        cazy = annotation_wrapper.CazyAnalysis(cazy_output=cazy_output, sql_path=sql_path)
        cazy.run()
    except IOError:
        logger.debug("NO cazy_output file!")


def insert_merops(merops_output, summary_sql_path):
    try:
        merops = annotation_wrapper.MeropsAnalysis(merops_output=merops_output, summary_sql_path=summary_sql_path)
        merops.run()
    except IOError:
        logger.debug("No merops_output file!")


def insert_query_length(task_name, summary_sql_path, fasta_sql_path):
    try:
        query_length = annotation_wrapper.QueryAnalysis(task_name=task_name, summary_sql_path=summary_sql_path,
                                                        fasta_sql_path=fasta_sql_path
                                                        )
        query_length.run()
    except Exception as e:
        logger.error("Got an exception!", e)


def process_sequence():
    work_dir = args.work_dir
    task_name = args.task_name

    task_dir = os.path.join(work_dir, task_name)
    if not os.path.exists(task_dir):
        logger.info('Create dir', task_dir)
        os.mkdir(task_dir)
    summary_sql_path = os.path.join(task_dir, 'qprotein_results.db')
    uniprot_db = os.path.join(work_dir, 'qprotein_db.db')
    sprot_dmnd_output = os.path.join(task_dir, 'sprot_70_200_90.out')
    trembl_dmnd_output = os.path.join(task_dir, 'trembl_70_200_90.out')
    cazy_output = os.path.join(task_dir, 'cazy_overview.txt')
    merops_output = os.path.join(task_dir, 'merops.out')

    create_sql(summary_sql_path=summary_sql_path)
    insert_sprot(dmnd_output=sprot_dmnd_output, summary_sql_path=summary_sql_path)
    annotate_sprot(uniprot_db=uniprot_db, summary_sql_path=summary_sql_path)
    insert_trembl(dmnd_output=trembl_dmnd_output, summary_sql_path=summary_sql_path)
    annotate_trembl(summary_sql_path=summary_sql_path, uniprot_db=uniprot_db)
    insert_cazy(cazy_output=cazy_output, sql_path=summary_sql_path)
    insert_merops(merops_output=merops_output, summary_sql_path=summary_sql_path)
    insert_query_length(task_name=task_name, summary_sql_path=summary_sql_path, fasta_sql_path=uniprot_db)


if __name__ == '__main__':
    process_sequence()
