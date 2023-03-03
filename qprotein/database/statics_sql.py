# -*- coding: utf-8 -*-
# @Time    : 2022/12/1 19:06
# @Author  : Zhixin Dou

import sqlite3
import os

def count_func(sql_db):
    connect = sqlite3.connect(sql_db)
    cur = connect.cursor()
    sql = "SELECT UniDesc, COUNT (UniDesc) FROM AFmeta GROUP BY UniDesc ORDER BY COUNT (UniDesc) DESC LIMIT 26"
    ret = cur.execute(sql)
    results = ret.fetchall()
    for i in results:
        print(i)


    with open(r'D:\subject\active\PyMulstruct\data\results_for_function_cluster.txt', 'w') as wt:
        for i in results:
            sql_output = "SELECT * FROM AFmeta WHERE UniDesc=? Limit 1"
            cur.execute(sql_output, (i[0],))
            ret_output = cur.fetchall()
            for result in ret_output:
                str_result = [str(i) for i in result]
                wt.write('\t'.join(str_result)+'\n')
    cur.close()
    connect.close()

sql_db = r'G:\DB\AFmeta.db'
count_func(sql_db=sql_db)
