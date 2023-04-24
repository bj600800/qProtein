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
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)

parser = argparse.ArgumentParser(description='Crawl query structure with multiprocessing.')
parser.add_argument('--work_dir', required=True, help='Specify the work directory path')
parser.add_argument('--db_dir', required=True, help='Specify all of the database directory path')
args = parser.parse_args()


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
    db_sql_path = os.path.join(work_dir, 'qprotein_db.db')
    sprot_dat_file = os.path.join(db_dir, "uniprot", "annotation", "uniprot_sprot.dat.gz")
    trembl_dat_file = os.path.join(db_dir, "uniprot", "annotation", "uniprot_trembl.dat.gz")

    create_sprot_dat_table(db_sql_path=db_sql_path, sprot_dat_file=sprot_dat_file)
    create_trembl_dat_table(db_sql_path=db_sql_path, trembl_dat_file=trembl_dat_file)
    logger.info(f"Path to SQL database: {os.path.abspath(db_sql_path)}")


run()
