# -*- coding: utf-8 -*-
# @Time    : 2023/2/24 11:34
# @Author  : Zhixin Dou
"""
Base class for SQL builder class
"""
from abc import abstractmethod
import sqlite3


class SqlBuilder(object):

    def __init__(self, sql_db, table_name, column_definition):
        self.sql_db = sql_db
        self.table_name = table_name
        self.column_definition = column_definition
        self.records = []

    @abstractmethod
    def read_file(self):
        """
        This is an abstract method that needs to be implemented in the subclass.
        """
        raise NotImplementedError("The read method must be implemented by a subclass")

    def create_cath_sql(self):
        connect = sqlite3.connect(self.sql_db)
        cursor = connect.cursor()
        cursor.execute(f"DROP TABLE IF EXISTS {self.table_name}")
        definition_str = ','.join([' '.join(i) for i in self.column_definition.items()])
        cursor.execute(f"CREATE TABLE {self.table_name} ({definition_str})")
        return cursor

    @abstractmethod
    def pack_records(self, desc_content):
        """
        This is an abstract method that needs to be implemented in the subclass.
        """
        raise NotImplementedError("The read method must be implemented by a subclass")

    def insert_batch(self, cursor, records: list):
        cursor.execute("begin")
        cursor.executemany(f"INSERT INTO {self.table_name} VALUES(?, ?, ?, ?, ?, ?)", records)
        cursor.execute("commit")

    def create_index(self, cursor):
        column_name = ', '.join(self.column_definition.keys())
        cursor.execute(f"CREATE INDEX cath_desc_idx ON cath_desc ({column_name})")

    @abstractmethod
    def run(self):
        """
        This is an abstract method that needs to be implemented in the subclass.
        """
        raise NotImplementedError("The read method must be implemented by a subclass")
