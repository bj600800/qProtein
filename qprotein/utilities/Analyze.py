"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/04/16

# Description: data statistics
# ------------------------------------------------------------------------------
"""
import os
import re

from tqdm import tqdm
from qprotein.database.sqlite3_builder import SqlBuilder
from qprotein.database.sqlite3_searcher import SqlSearch

class Statistics:
    def __init__(self, high_dir, good_dir, summary_sql, topt_results, task_name, plddt_output, topt_output, statistic_output):
        self.high_dir = high_dir
        self.good_dir = good_dir
        self.summary_sql = summary_sql
        self.topt_results = topt_results
        self.task_name = task_name
        self.plddt_output = plddt_output
        self.topt_output = topt_output
        self.statistic_output = statistic_output

    def get_file(self, cif_dir):
        files = os.listdir(cif_dir)
        return files

    def get_query_name(self, file):
        query_name = os.path.splitext(file)[0].split("#")[0]
        return query_name

    def get_plddt(self, file):
        plddt = os.path.splitext(file)[0].split("#")[2]
        return plddt

    def get_topt(self):
        with open(self.topt_results, "r") as rf:
            topt = [line.split("\t")[1:] for line in rf.readlines()[1:] if line.split("\t")[1] != "name"]
        print(topt)
        return topt

    def count_cif(self):
        count = len(os.listdir(self.high_dir)) + len(os.listdir(self.good_dir))
        return count

    def get_high_prop(self):
        count = self.count_cif()
        high_prop = len(os.listdir(self.high_dir))/count
        return high_prop

    def get_good_prop(self):
        count = self.count_cif()
        good_prop = len(os.listdir(self.good_dir))/count
        return good_prop


class SearchParam(SqlBuilder, SqlSearch):
    # TODO 结果分析应该从results_summary表中查询，而不是重新构建一个表。
    # 出现重复是因为之前 Tibet[:300000], [300000:]出现了重叠。重叠部分将近8万个结构。
    def __init__(self, param, task_dir):
        super(SearchParam, self).__init__(sql_db=os.path.join(task_dir, "structure_results.db"),
                                          table_name="results",
                                          column_definition=[("query_name", "TEXT"), ("helix", "FLOAT"), ("sheet", "FLOAT"),
                                                             ("loop", "FLOAT"), ("turn", "FLOAT"),
                                                             ("bend", "FLOAT"), ("bridge", "FLOAT"),
                                                             ("hbond_density", "FLOAT"), ("hbond_avg_energy", "FLOAT"),
                                                             ("apolar", "FLOAT"), ("polar", "FLOAT"),
                                                             ("positive", "FLOAT"), ("negative", "FLOAT"),
                                                             ("polar_area", "FLOAT"), ("apolar_area", "FLOAT")])
        self.sql_db = os.path.join(task_dir, "structure_results.db")
        self.param = param  # 70_200_90
        self.task_dir = task_dir
        self.structure_results_file = os.path.join(task_dir, "structure_results.csv")
        self.structure_dir = os.path.join(task_dir, "high_plddt")
        self.target_structure_dir = os.path.join(task_dir, "uniprot_grid_output", param)
        if not os.path.exists(self.target_structure_dir):
            os.mkdir(self.target_structure_dir)

    def read_out(self):
        uniprot_out_files = ["sprot_"+self.param+".out", "trembl_"+self.param+".out"]
        files = [os.path.join(self.task_dir, "uniprot_grid_output", i) for i in uniprot_out_files]
        query_name = []
        for file in files:
            with open(file, "r") as f:
                content = f.readlines()
            for line in content:
                line = line.split("\t")
                query_name.append(line[0])
        return query_name

    def build_sql_results(self):
        sql = f"CREATE TABLE {self.table_name} (query_name TEXT)"
        cursor = self.create_table(table_name=self.table_name, sql_db=self.sql_db, sql=sql)
        self.add_column(cursor=cursor, table_name=self.table_name, column_definition=self.column_definition)
        columns = [i[0] for i in self.column_definition]
        with open(self.structure_results_file, "r") as f:
            records = [(i.rstrip().split(",")[:-3]) for i in f.readlines()][1:]
        self.insert_many_columns(cursor=cursor, columns=columns, records=records, values_num=15)
        self.create_index(cursor=cursor, columns=columns)

    def fetch_sql_data(self, query_name):
        cursor = self.connect_sql(sql_db=self.sql_db)
        batch_size = 100000
        data = []
        total_query = len(query_name)
        for i in range(0, total_query, batch_size):
            batch_query = query_name[i:i + batch_size]
            sql_cmd = "SELECT * FROM {} WHERE query_name IN ({})" \
                      .format(self.table_name, ', '.join(['?'] * len(query_name)))
            results_list = self.fetch_many_results(cursor=cursor, sql_cmd=sql_cmd, search_targets=batch_query)
            data.extend(results_list)
        return data

    def get_structure_attr(self):
        if not os.path.exists(self.sql_db):
            self.build_sql_results()
        query_name = self.read_out()
        data = self.fetch_sql_data(query_name=query_name)
        print(data)
        input()
        return data

