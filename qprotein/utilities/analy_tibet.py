# -*- coding: utf-8 -*-
# @Time    : 2022/12/2 12:53
# @Author  : Zhixin Dou
import os
from collections import defaultdict, Counter
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def analy(txt):
    sample_dict = defaultdict(list)
    output_dict = {}
    with open(txt, 'r') as rd:
        content = [i.strip('\n').split() for i in rd.readlines()]
        for i in content:
            place = i[0].split('_')[0]
            sample_dict[place].append(i[1])
    for k, v in sample_dict.items():
        value_dict = {}
        for key, val in Counter(v).items():
            value_dict[key] = val
        output_dict[k] = value_dict

    df = pd.DataFrame(output_dict).fillna(value=0)
    print(df)
    df.to_csv(r'D:\subject\active\PyMulstruct\data\tibet.csv', sep=',')
    #         value_dict[i[1]].append(i[0])
    #     output_dict[k] = value_dict
    # print(output_dict['A2'])
    sns.heatmap(df, annot=False)
    plt.show()

txt = r'D:\subject\active\PyMulstruct\data\tibet.txt'
analy(txt)