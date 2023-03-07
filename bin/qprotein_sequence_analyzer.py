"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: perform sequence annotation and sql builder
# ------------------------------------------------------------------------------
"""
from qprotein.database import sqlite3_builder
from qprotein.database import sqlite3_searcher


class SeqAnalysis(sqlite3_builder, sqlite3_searcher):
    def __init__(self, sql_db):
        super.__init__(sql_db)

    def create_results_sql(self):
