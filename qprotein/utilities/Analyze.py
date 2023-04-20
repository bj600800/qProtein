"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/04/16

# Description: data statistics
# ------------------------------------------------------------------------------
"""
import os
from tqdm import tqdm


class Analyze:
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

    def run(self):
        files = self.get_file(self.high_dir)
        # with open(self.topt_output, 'w') as wf_t:
        #     topt = self.get_topt()
        #     for t in topt:
        #         wf_t.write(t[1])

        with open(self.plddt_output, "w") as wf_p:
            for file_name in files:
                plddt = self.get_plddt(file_name)
                wf_p.write(plddt+"\n")

        with open(self.statistic_output, "w") as wf_s:
            high_prop = self.get_high_prop()
            good_prop = self.get_good_prop()
            wf_s.write("high_prop:\t"+str(high_prop)+"\n")
            wf_s.write("good_prop:\t"+str(good_prop)+"\n")


if __name__ == '__main__':
    task_name = "manure"
    root_path = r"D:\subject\active\1-qProtein\data"
    work_path = os.path.join(root_path, task_name)
    summary_sql = os.path.join(work_path, "qprotein_results.db")
    topt_results = os.path.join(work_path, task_name+"_high_plddt.tsv")
    high_dir = os.path.join(work_path, 'high_plddt')
    good_dir = os.path.join(work_path, 'good_plddt')
    plddt_output = os.path.join(work_path, "plddt_output.txt")
    topt_output = os.path.join(work_path, "topt_output.txt")
    statistic_output = os.path.join(work_path, "statistic.txt")
    task = Analyze(high_dir=high_dir, good_dir=good_dir, summary_sql=summary_sql,
                   topt_output=topt_results, task_name=task_name, plddt_output=plddt_output,
                   topt_results=topt_output, statistic_output=statistic_output)
    task.run()
