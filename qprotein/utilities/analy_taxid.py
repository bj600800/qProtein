# -*- coding: utf-8 -*-
# @Time    : 2022/11/28 9:38
# @Author  : Zhixin Dou
"""
Convert file format
From:
Organism name; prediction number; very_low_prop; low_

"""
import os
import csv
from collections import defaultdict
import numpy as np
import pandas as pd


def read_csv(csv_file):
    phly_dic = defaultdict(list)
    csv_reader = csv.reader(open(csv_file))
    for line in csv_reader:
        list_per_line = []
        list_per_line.append(line[6])
        list_per_line.append(line[5])
        list_per_line.append(line[4])
        list_per_line.append(line[3])
        list_per_line.append(line[2])
        phly_dic[line[8].split('__')[1]].append(list_per_line)
    return phly_dic

def calc_dic(dic):
    out_dic = defaultdict(dict)
    for phly, num_plddt in dic.items():
        num = []
        high = []
        confi = []
        low = []
        very_low = []

        for i in num_plddt:
            num.append(int(i[0]))
            high.append(float(i[1]))
            confi.append(float(i[2]))
            low.append(float(i[3]))
            very_low.append(float(i[4]))

        # print(phly)
        # print('number:', sum(np.array(num)))
        # print('high:', np.array(high))
        # print('average high:', np.average(high))
        # print('weights_average:', int(np.average(num, weights=high)))
        # input()
        out_dic[phly]['Total_predictions'] = int(sum(num))
        out_dic[phly]['high_prop'] = np.average(np.array(high))
        out_dic[phly]['confi_prop'] = np.average(np.array(confi))
        out_dic[phly]['low_prop'] = np.average(np.array(low))
        out_dic[phly]['very_low_prop'] = np.average(np.array(very_low))
    df = pd.DataFrame(out_dic).T
    print(df)
    return df

def main():
    file_dir = r'D:\subject\active\PyMulstruct\alphafold\alphafold_metadata\4-class\complete'
    file_list = os.listdir(file_dir)
    output_dir = r'D:\subject\active\PyMulstruct\alphafold\alphafold_metadata\output'
    for i in file_list:
        csv_file = os.path.join(file_dir, i)
        phly_dic = read_csv(csv_file)
        out_dic = calc_dic(phly_dic)
        output_csv_path = os.path.join(output_dir, '{i}_statistic.csv'.format(i=i.split('.')[0]))
        out_dic.to_csv(output_csv_path)

main()
