import pandas as pd
import numpy as np
import os
from collections import defaultdict

phylum = []
bio_class = []
order = []
family = []
genus = []
species = []


def append_lineage(frame, txt):
    frame = frame.values.tolist()
    for line_meta in frame:
        n = 0
        for line_txt in txt.values.tolist():
            if n < 1:  # multi to one
                if line_meta[1] == int(line_txt[0]):
                    n += 1
                    # line_meta.append(line_txt[1])
                    for i in line_txt[1:]:
                        line_meta.append(i)

    frame = pd.DataFrame(frame)
    frame.dropna(axis=1, how='all')
    return frame


def df2excel(excel_file, text_file, rd_sheet_name, wt_sheet_name):
    frame = pd.read_excel(excel_file, sheet_name=rd_sheet_name)
    txt = pd.read_table(text_file, sep='\t')
    frame = append_lineage(frame, txt)
    columns_name = 'name taxID mean_plddt_very_low mean_plddt_low mean_plddt_confident mean_plddt_very_high num_predictions lineage'.split(
        ' ')
    print('Writing to csv')
    data = pd.DataFrame(frame)
    # pd.set_option('display.max_rows', None) # show all columns
    print(data)
    data.to_csv(wt_sheet_name+'.csv', index=False, header=False)


def main():
    taxid_root_path = r'D:\subject\active\PyMulstruct\alphafold\alphafold_metadata\4-class'
    meta_xls_path = r'D:\subject\active\PyMulstruct\alphafold\alphafold_metadata\groupByorganism-taxid.xlsx'
    for i in range(1, 2):
        print('processing file:', 'lineage_class{i}.txt'.format(i=i))
        txt_file = os.path.join(taxid_root_path, 'lineage_class{num}.txt'.format(num=i))
        rd_sheet_name = str(i)
        wt_sheet_name = os.path.join(taxid_root_path, 'lineage_class{num}'.format(num=i))
        df2excel(meta_xls_path, txt_file, rd_sheet_name, wt_sheet_name)


main()
