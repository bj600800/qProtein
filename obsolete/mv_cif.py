"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/15

# Description: 
# ------------------------------------------------------------------------------
"""
import os
import sqlite3
import shutil


def connect_sql(sql_db):
    connect = sqlite3.connect(sql_db)
    cursor = connect.cursor()
    return cursor


def target_cif(dir_path):
    file_list = os.listdir(dir_path)
    return file_list


def get_cazy_family(cursor, file_list):
    sql = 'SELECT query_name, cazy_family FROM results_summary WHERE query_name in ({})'\
        .format(', '.join(['?'] * len(file_list)))
    query_list = ['_'.join(file.split('.')[0].split('_')[:6]) for file in file_list]
    cursor.execute(sql, query_list)
    return_dat = cursor.fetchall()
    return return_dat


def make_dir(new_dir):
    if not os.path.exists(new_dir):
        print(f'Making directory in {new_dir}')
        os.makedirs(os.path.join(new_dir))


def move_cif(return_dat, dir_path, new_dir):
    make_dir(new_dir)
    file_list = os.listdir(dir_path)
    for file in file_list:
        for dat in return_dat:
            if dat[0] == '_'.join(file.split('.')[0].split('_')[:6]):
                src_file = os.path.join(dir_path, file)
                src_file = open(src_file, 'rb')
                family = dat[1].split('(')[0]
                family_dir = os.path.join(new_dir, family)
                make_dir(family_dir)
                tgt_file = os.path.join(family_dir, file)
                tgt_file = open(tgt_file, 'wb')
                shutil.copyfileobj(src_file, tgt_file)



if __name__ == '__main__':
    sql_db = r'D:\subject\active\1-qProtein\data\tibet\qprotein_results.db'
    dir_path = r'D:\subject\active\1-qProtein\data\tibet\ident90_segment'
    new_dir = r'D:\subject\active\1-qProtein\data\tibet\ident90_new'
    cursor = connect_sql(sql_db)
    file_list = target_cif(dir_path)
    return_dat = get_cazy_family(cursor, file_list)
    move_cif(return_dat, dir_path, new_dir)
