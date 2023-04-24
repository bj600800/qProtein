"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/04/24

# Description: Add query fasta to qprotein_db
# ------------------------------------------------------------------------------
"""
import os.path
import sys
import argparse
from qprotein.database.query_seq_builder import FastaSql
from qprotein.utilities import logger

sys.path.append("..")
logger = logger.setup_log(name=__name__)

parser = argparse.ArgumentParser(description='Add query fasta sequences to qprotein_db')
parser.add_argument('--work_dir', required=True, help='Specify the work directory path')
parser.add_argument('--query_fasta', required=True, help='Path to the query fasta file')
args = parser.parse_args()


def create_query_seq_table(db_sql_path, fasta_file, task_name):
    try:
        fasta_db = FastaSql(sql_path=db_sql_path, fasta_file=fasta_file, table_name=task_name)
        fasta_db.run()
    except Exception as e:
        logger.debug(e)
        raise


def run():
    work_dir = args.work_dir
    fasta_file = args.query_fasta
    task_name = os.path.splitext(os.path.basename(fasta_file))[0]
    db_sql_path = os.path.join(work_dir, 'qprotein_db.db')
    create_query_seq_table(db_sql_path=db_sql_path, fasta_file=fasta_file, task_name=task_name)
    logger.info(f"Add query fasta to SQL database: {os.path.abspath(db_sql_path)}")


run()

