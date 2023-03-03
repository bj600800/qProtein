# -*- coding: utf-8 -*-
# @Time    : 2022/12/19 14:32
# @Author  : Zhixin Dou
"""
sqlite3 search targets against uniprot_db
"""
import sqlite3
import time


def search_sql(uniprot_sql_db, query_file, output_file):
    start_time = time.time()
    connect = sqlite3.connect(uniprot_sql_db)
    cur = connect.cursor()
    with open(output_file, 'a+') as wt:
        with open(query_file, 'r') as rd:
            target_list = [i.split(' ')[0] for i in rd.readlines()]
        load_time = time.time()
        print('Targets loaded time consuming: {:.2f} seconds'.format(load_time - start_time))
        print('-' * 30)
        counter = 0
        name_list = []
        for i in target_list:
            if counter == 10000:
                print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
                search_begin_time = time.time()
                sql = "SELECT * FROM uniprot INDEXED BY idx WHERE NAME in ({name})".format(
                    name=','.join(['?'] * len(name_list)))
                cur.execute(sql, name_list)
                row = cur.fetchall()
                count = len(row)
                if row:
                    for fasta in row:
                        wt.write('>' + fasta[0] + '\n')
                        wt.write(fasta[1] + '\n')
                        wt.flush()
                    print('Write target fasta sequences', count)
                    end_time = time.time()
                    print("Time consuming for 10000 sequences searching: {:.2f} seconds".format(
                        end_time - search_begin_time))
                    print('-' * 30)
                counter = 0
                name_list = []
            else:
                name_list.append(i)
                counter += 1

        if counter > 0:
            print('process the remains.')
            sql = "SELECT * FROM uniprot WHERE NAME in ({name})".format(
                name=','.join(['?'] * len(name_list)))
            cur.execute(sql, name_list)
            row = cur.fetchall()
            if row:
                with open(output_file, 'a+') as wt_remain:
                    for fasta in row:
                        wt_remain.write('>' + fasta[0] + '\n')
                        wt_remain.write(fasta[1] + '\n')
                        wt_remain.flush()
                    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))

    print('All sequences done!')
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))


uniprot_sql_db = r'G:/DB/uniprot_fasta.db'
query_file = r'G:/DB/continue.txt'
output_file = r'G:/DB/alphafold_mean90_seqs.fasta'

search_sql(uniprot_sql_db, query_file, output_file)
