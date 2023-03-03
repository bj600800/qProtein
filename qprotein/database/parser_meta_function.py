# -*- coding: utf-8 -*-
# @Time    : 2022/11/30 16:05
# @Author  : Zhixin Dou

import csv
import os
import sqlite3


def create_AFmeta_sql(AFmeta_sql_db, meta_csv_list):
    if os.path.exists(AFmeta_sql_db):
        return
    connect = sqlite3.connect(AFmeta_sql_db)
    cur = connect.cursor()
    cur.execute("DROP TABLE IF EXISTS AFmeta")
    cur.execute('''CREATE TABLE AFmeta \
    (AFID TEXT, UniACC TEXT, UniDesc TEXT, Reviewed TEXT, StartRes INTEGER, EndRes INTEGER, mean_plddt REAL,\
    FracVeryLow REAL, FracLow REAL, FracConfi REAL, FracVeryHigh REAL)''')
    counter = 0
    record = []
    for meta_csv in meta_csv_list:
        csv_reader = csv.reader(open(meta_csv))
        for line in csv_reader:
            if counter == 500000:
                cur.execute('begin')
                cur.executemany('INSERT INTO AFmeta VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', record)
                cur.execute('commit')
                print('Write fasta files', counter)
                counter = 0
                record = []
            else:
                record.append(tuple(line))
                counter += 1


def search_sql(uniprot_sql_db, target_file, output):
    connect = sqlite3.connect(uniprot_sql_db)
    cur = connect.cursor()
    with open(output, 'w') as wt:
        with open(target_file, 'r') as rd:
            target_list = [i[0] for i in rd.readlines()]
        for i in target_list:
            cur.execute("SELECT * FROM uniprot WHERE name=?", (i,))
            row = cur.fetchone()
            wt.write('>' + row[0] + '\n')
            wt.write(row[1] + '\n')


AFmeta_sql_db = r'G:\DB\AFmeta.db'
meta_root = r'G:\DB\alphafold_plddt_mean_greater90_meta'
meta_csv_list = [os.path.join(meta_root, i) for i in os.listdir(meta_root)]
create_AFmeta_sql(AFmeta_sql_db, meta_csv_list)
# search_sql(uniprot_sql_db, txt_file, output)
