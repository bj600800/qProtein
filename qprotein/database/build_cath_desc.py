# -*- coding: utf-8 -*-
# @Time    : 2023/2/24 11:34
# @Author  : Zhixin Dou
"""
Build cath description sqlite file
"""

import sqlite3
import time


def read_desc_file(desc_file):
    with open(desc_file, 'r') as f:
        desc = [line for line in f.readlines() if not line.startswith('#')]
    return desc


def create_cath_sql(sql_db):
    connect = sqlite3.connect(sql_db)
    cursor = connect.cursor()
    cursor.execute("DROP TABLE IF EXISTS cath_desc")
    cursor.execute("CREATE TABLE cath_desc "
                   "(domain_name TEXT, cath_code TEXT, cath_class TEXT, "
                   "cath_architecture TEXT, cath_topology TEXT, cath_homologous TEXT)")


def build_sqlite(desc, sql_db):
    connect = sqlite3.connect(sql_db)
    cursor = connect.cursor()

    domain_name = ''
    cath_code = ''
    cath_class = ''
    cath_architecture = ''
    cath_topology = ''
    cath_homologous = ''

    records = []
    record_counter = 0
    for line in desc:
        _list = []
        if line.startswith("DOMAIN"):
            # GET accession number, the first one
            domain_name = line.split()[1]

        elif line.startswith("CATHCODE"):
            cath_code = line.split()[1]

        elif line.startswith("CLASS"):
            cath_class = line.split()[1]

        elif line.startswith("ARCH"):
            cath_architecture = line.split()[1]

        elif line.startswith("TOPOL"):
            cath_topology = line.split()[1]

        elif line.startswith("HOMOL"):
            cath_homologous = line.split()[1]

        elif "//\n" in line:
            records.append(tuple([domain_name, cath_code, cath_class, cath_architecture,
                                  cath_topology, cath_homologous]))
            record_counter += 1
            domain_name = ''
            cath_code = ''
            cath_class = ''
            cath_architecture = ''
            cath_topology = ''
            cath_homologous = ''

        if record_counter == 50000:
            cursor.execute("begin")
            cursor.executemany("INSERT INTO cath_desc VALUES(?, ?, ?, ?, ?, ?)", records)
            cursor.execute("commit")
            print('Insert 50000 items.')
            records = []
            record_counter = 0

    if record_counter > 0:
        print('Procession remaining records.')
        cursor.execute("begin")
        cursor.executemany(f"INSERT INTO cath_desc VALUES(?, ?, ?, ?, ?, ?)", records)
        cursor.execute("commit")
    print(f'Creating index for cath_desc.')
    cursor.execute(f"CREATE INDEX cath_desc_idx ON cath_desc (domain_name, cath_code, cath_class,"
                   "cath_architecture, cath_topology, cath_homologous)")

    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Finished sql_db building for cath_desc.')


def main():
    sql_db = r'G:\DB\cath\cath-classification-data_txt\cath.db'
    desc_file = r'G:\DB\cath\cath-classification-data_txt\cath-domain-description-file.txt'
    desc = read_desc_file(desc_file)
    create_cath_sql(sql_db)
    build_sqlite(desc, sql_db)


if __name__ == '__main__':
    main()
