"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: Service for building sql database
# ------------------------------------------------------------------------------
"""
import os.path
import sys
import argparse
sys.path.append("..")

from qprotein.database.uniprot_dat_builder import UniprotSql
from qprotein.database.query_seq_builder import FastaSql
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)

parser = argparse.ArgumentParser(description='Crawl query structure with multiprocessing.')
parser.add_argument('--work_dir', required=True, help='Specify the work directory path')
parser.add_argument('--db_dir', required=True, help='Specify all of the database directory path')
parser.add_argument('--only_query', required=False, help='True or False. Only create query sequence table')
parser.add_argument('--fasta_file', required=True, help='Path to the query fasta file')
args = parser.parse_args()


def create_query_seq_table(db_sql_path, fasta_file, task_name):
    try:
        fasta_db = FastaSql(sql_path=db_sql_path, fasta_file=fasta_file, table_name=task_name)
        fasta_db.run()
    except Exception as e:
        logger.debug(e)
        raise


def create_sprot_dat_table(db_sql_path, sprot_dat_file):
    try:
        sprot_builder = UniprotSql(dat_file=sprot_dat_file, sql_path=db_sql_path, table_name='sprot_dat')
        sprot_builder.run()
    except Exception as e:
        logger.debug(e)
        raise

def create_trembl_dat_table(db_sql_path, trembl_dat_file):
    try:
        sprot_builder = UniprotSql(dat_file=trembl_dat_file, sql_path=db_sql_path, table_name='trembl_dat')
        sprot_builder.run()
    except Exception as e:
        logger.debug(e)
        raise

def run():
    work_dir = args.work_dir
    db_dir = args.db_dir
    only_query = args.only_query
    fasta_file = args.fasta_file

    task_name = os.path.splitext(os.path.basename(fasta_file))[0]
    db_sql_path = os.path.join(work_dir, 'qprotein_db.db')
    # sprot_dat_file = os.path.join(db_dir, "uniprot", "annotation", "uniprot_sprot.dat.gz")
    sprot_dat_file = os.path.join(db_dir, "uniprot", "annotation", "test.dat.gz")
    trembl_dat_file = os.path.join(db_dir, "uniprot", "annotation", "uniprot_trembl.dat.gz")

    if only_query:
        create_query_seq_table(db_sql_path=db_sql_path, fasta_file=fasta_file, task_name=task_name)
        return
    else:
        create_query_seq_table(db_sql_path=db_sql_path, fasta_file=fasta_file, task_name=task_name)
        create_sprot_dat_table(db_sql_path=db_sql_path, sprot_dat_file=sprot_dat_file)
        #create_trembl_dat_table(db_sql_path=db_sql_path, trembl_dat_file=trembl_dat_file)

    logger.info(f"Path to SQL database: {os.path.abspath(db_sql_path)}")

run()
