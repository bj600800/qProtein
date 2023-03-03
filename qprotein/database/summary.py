# -*- coding: utf-8 -*-
# @Time    : 2023/1/3 14:35
# @Author  : Zhixin Dou
"""
parse and summary all annotation outputs into sqlite database.
"""
import os
import sqlite3
import time


def create_sqlite(sql_db):
    if os.path.exists(sql_db):
        print(sql_db, ' existed, deleting...')
        os.remove(sql_db)
    connect = sqlite3.connect(sql_db)
    cursor = connect.cursor()
    cursor.execute("DROP TABLE IF EXISTS summary")
    cursor.execute("CREATE TABLE summary "
                   "(acc_id TEXT, ec_number, cazy_family TEXT,"
                   "merops70_family TEXT, sprot100_acc TEXT, sprot_name TEXT,"
                   "sprot_ec_number TEXT, sprot_go_component TEXT, sprot_go_process TEXT,"
                   "sprot_go_function TEXT, sprot_interpro TEXT, sprot_pfam TEXT,"
                   "sprot_start INTEGER, sprot_end INTEGER,"
                   "trembl100_acc TEXT, trembl_name TEXT, "
                   "trembl_ec_number TEXT, trembl_go_component TEXT, trembl_go_process TEXT,"
                   "trembl_go_function TEXT, trembl_interpro TEXT, trembl_pfam TEXT,"
                   "trembl_start INTEGER, trembl_end INTEGER, "
                   "alphafold_mmcif_file TEXT, PDB100_acc TEXT, pdb_mmcif_file TEXT)")


def parse_cazy(cazy_output, sql_db):
    """
    priority: 1st, 2nd, 3rd columns
    """
    time1 = time.time()
    print('=' * 30)
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Starting to parse CAZy output, please wait...')
    connect = sqlite3.connect(sql_db)
    cursor = connect.cursor()
    with open(cazy_output, 'r') as rd:
        text = [i.split('\t') for i in rd.readlines()[1:]]
    records = []
    record_counter = 0
    for i in text:
        if record_counter == 10000:
            cursor.execute("begin")
            cursor.executemany("INSERT INTO summary(acc_id, ec_number, cazy_family) VALUES(?, ?, ?)", records)
            cursor.execute("commit")
            records = []
            record_counter = 0

        _dict = {}
        _dict['acc_id'] = i[0]
        if i[1] == '-':
            _dict['ec_number'] = ''
        else:
            _dict['ec_number'] = i[1]
        if i[2] != '-':
            _dict['cazy_family'] = i[2].split('(')[0]
        elif i[3] != '-':
            _dict['cazy_family'] = i[3]
        elif i[4] != '-':
            _dict['cazy_family'] = i[4]
        records.append(tuple(
            [_dict['acc_id'], _dict['ec_number'], _dict['cazy_family']]))
        record_counter += 1

    if record_counter > 0:
        cursor.execute("begin")
        cursor.executemany("INSERT INTO summary(acc_id, ec_number, cazy_family) VALUES(?, ?, ?)", records)
        cursor.execute("commit")
    cursor.execute("CREATE INDEX summary_idx_cazy ON summary (acc_id, ec_number, cazy_family)")

    time2 = time.time()
    time_consum = time2 - time1
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Finished parsing cazy output in {:.2f} secs.'.format(time_consum))


def parse_merops(merops_output, sql_db):
    time1 = time.time()
    print('=' * 30)
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Starting to parse merops output, please wait...')
    connect = sqlite3.connect(sql_db)
    cursor = connect.cursor()
    with open(merops_output, 'r') as rd:
        text = [i.split('\t') for i in rd.readlines()]
    record_counter = 0
    records = []
    for i in text:
        if record_counter == 10000:
            cursor.execute("begin")
            cursor.executemany("REPLACE INTO summary(acc_id, merops70_family) VALUES(?,?)", records)
            cursor.execute("commit")
            record_counter = 0
            records = []

        # For the N+1 item, appending to list as well.
        _dict = {}
        _dict['acc_id'] = i[0]
        _dict['merops70_family'] = i[1].split('#')[1]
        records.append(tuple([_dict['acc_id'], _dict['merops70_family']]))
        record_counter += 1

    if record_counter > 0:
        print('Processing remaining tasks.')
        cursor.execute("begin")
        cursor.executemany("REPLACE INTO summary(acc_id, merops70_family) VALUES(?,?)", records)
        cursor.execute("commit")
    cursor.execute("CREATE INDEX summary_idx_merops ON summary (merops70_family)")

    time2 = time.time()
    time_consum = time2 - time1
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Finished parsing cazy output in {:.2f} secs.'.format(time_consum))


def parse_sprot_dmnd(sprot_output, sql_db):
    time1 = time.time()
    print('=' * 30)
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Starting to parse Swiss-prot output, please wait...')
    connect = sqlite3.connect(sql_db)
    cursor = connect.cursor()
    with open(sprot_output, 'r') as rd:
        text = [i.split('\t') for i in rd.readlines()]
    record_counter = 0
    records = []
    for i in text:
        if record_counter == 10000:
            cursor.execute("begin")
            cursor.executemany("UPDATE summary SET sprot100_acc = ?, sprot_start = ?, sprot_end = ?, where acc_id=?",
                               records)
            cursor.execute("commit")
            record_counter = 0
            records = []

        _dict = {}
        _dict['acc_id'] = i[0]
        _dict['sprot100_acc'] = i[1].split('|')[1]
        _dict['sprot_start'] = int(i[8])
        _dict['sprot_end'] = int(i[9])
        records.append(tuple([_dict['sprot100_acc'], _dict['sprot_start'], _dict['sprot_end'], _dict['acc_id']]))
        record_counter += 1

    if record_counter > 0:
        print('Processing remaining tasks.')
        cursor.execute("begin")
        cursor.executemany("UPDATE summary SET sprot100_acc = ?, sprot_start = ?, sprot_end = ? where acc_id=?",
                           records)
        cursor.execute("commit")
    cursor.execute("CREATE INDEX summary_idx_sprot ON summary (sprot100_acc)")

    time2 = time.time()
    time_consum = time2 - time1
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Finished parsing Swiss-prot output in {:.2f} secs.'.format(time_consum))


def annot_sprot_sql(sprot_db, sql_db):
    time1 = time.time()
    print('=' * 30)
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Starting to search sprot sqlite, please wait...')

    # fetch acc from sql_db in table (summary)
    connect_sql_db = sqlite3.connect(sql_db)
    cursor_sql_db = connect_sql_db.cursor()
    cursor_sql_db.execute("SELECT sprot100_acc FROM summary WHERE sprot100_acc != ''")
    ret_acc_list = [i[0] for i in cursor_sql_db.fetchall()]

    # fetch dat_info from sprot_db in table (uniprot_sprot)
    connect_sprot_db = sqlite3.connect(sprot_db)
    cursor_sprot_db = connect_sprot_db.cursor()
    cursor_sprot_db.execute("SELECT * FROM uniprot_sprot WHERE accession in ({acc_id})"
                            .format(acc_id=','.join(['?'] * len(ret_acc_list))), ret_acc_list)
    ret_dat_list = cursor_sprot_db.fetchall()

    # update
    record_counter = 0
    records = []
    for i in ret_dat_list:
        if record_counter == 10000:
            cursor_sql_db.execute("begin")
            cursor_sql_db.executemany(
                "UPDATE summary SET sprot_name = ?, sprot_ec_number = ?, sprot_go_component = ?,"
                "sprot_go_process = ?, sprot_go_function = ?, sprot_interpro = ?, sprot_pfam = ? where sprot100_acc=?",
                records)
            cursor_sql_db.execute("commit")
            record_counter = 0
            records = []

            # For the N+1 item, appending to list as well.
        _dict = {}
        _dict['sprot100_acc'] = i[0]
        _dict['sprot_name'] = i[1]
        _dict['sprot_ec_number'] = i[2]
        _dict['sprot_go_component'] = i[3]
        _dict['sprot_go_process'] = i[4]
        _dict['sprot_go_function'] = i[5]
        _dict['sprot_interpro'] = i[6]
        _dict['sprot_pfam'] = i[7]
        records.append(tuple([_dict['sprot_name'], _dict['sprot_ec_number'], _dict['sprot_go_component'],
                              _dict['sprot_go_process'], _dict['sprot_go_function'], _dict['sprot_interpro'],
                              _dict['sprot_pfam'], _dict['sprot100_acc']]))
        record_counter += 1

    if record_counter > 0:
        print('Processing remaining tasks.')
        cursor_sql_db.execute("begin")
        cursor_sql_db.executemany(
            "UPDATE summary SET sprot_name = ?, sprot_ec_number = ?, sprot_go_component = ?,"
            "sprot_go_process = ?, sprot_go_function = ?, sprot_interpro = ?, sprot_pfam = ? where sprot100_acc=?",
            records)
        cursor_sql_db.execute("commit")
    cursor_sql_db.execute("CREATE INDEX summary_idx_sprot_dat ON summary (sprot_ec_number, sprot_go_component,"
                          "sprot_go_process, sprot_go_function, sprot_interpro, sprot_pfam)")

    time2 = time.time()
    time_consum = time2 - time1
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Finished parsing sprot dat_file in {:.2f} secs.'.format(time_consum))


def parse_trembl_dmnd(trembl_output, sql_db):
    time1 = time.time()
    print('=' * 30)
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Starting to parse Trembl output, please wait...')
    connect = sqlite3.connect(sql_db)
    cursor = connect.cursor()
    with open(trembl_output, 'r') as rd:
        text = [i.split('\t') for i in rd.readlines()]
    record_counter = 0
    records = []
    for i in text:
        if record_counter == 10000:
            cursor.execute("begin")
            cursor.executemany("UPDATE summary SET trembl100_acc = ?, trembl_start = ?, trembl_end = ? where acc_id=?",
                               records)
            cursor.execute("commit")
            record_counter = 0
            records = []

            # For the N+1 item, appending to list as well.
        _dict = {}
        _dict['acc_id'] = i[0]
        _dict['trembl100_acc'] = i[1].split('|')[1]
        _dict['trembl_start'] = int(i[8])
        _dict['trembl_end'] = int(i[9])
        records.append(tuple([_dict['trembl100_acc'], _dict['trembl_start'], _dict['trembl_end'], _dict['acc_id']]))
        record_counter += 1

    if record_counter > 0:
        print('Processing remaining tasks.')
        cursor.execute("begin")
        cursor.executemany("UPDATE summary SET trembl100_acc = ?, trembl_start = ?, trembl_end = ? where acc_id=?",
                           records)
        cursor.execute("commit")
    cursor.execute("CREATE INDEX summary_idx_trembl ON summary (trembl100_acc)")

    time2 = time.time()
    time_consum = time2 - time1
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Finished parsing Trembl output in {:.2f} secs.'.format(time_consum))


def annot_trembl_sql(trembl_db, sql_db):
    time1 = time.time()
    print('=' * 30)
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Starting to search trembl sqlite, please wait...')

    # fetch acc from sql_db in table (summary)
    connect_sql_db = sqlite3.connect(sql_db)
    cursor_sql_db = connect_sql_db.cursor()
    cursor_sql_db.execute("SELECT trembl100_acc FROM summary WHERE trembl100_acc != ''")
    ret_acc_list = [i[0] for i in cursor_sql_db.fetchall()]

    # fetch dat_info from sprot_db in table (uniprot_sprot)
    connect_sprot_db = sqlite3.connect(trembl_db)
    cursor_sprot_db = connect_sprot_db.cursor()
    cursor_sprot_db.execute("SELECT * FROM uniprot_trembl WHERE accession in ({acc_id})"
                            .format(acc_id=','.join(['?'] * len(ret_acc_list))), ret_acc_list)
    ret_dat_list = cursor_sprot_db.fetchall()
    print(ret_dat_list)
    # update
    record_counter = 0
    records = []
    for i in ret_dat_list:
        if record_counter == 10000:
            cursor_sql_db.execute("begin")
            cursor_sql_db.executemany(
                "UPDATE summary SET trembl_name = ?, trembl_ec_number = ?, trembl_go_component = ?,"
                "trembl_go_process = ?, trembl_go_function = ?, trembl_interpro = ?, trembl_pfam = ? where trembl100_acc=?",
                records)
            cursor_sql_db.execute("commit")
            record_counter = 0
            records = []

            # For the N+1 item, appending to list as well.
        _dict = {}
        _dict['trembl100_acc'] = i[0]
        _dict['trembl_name'] = i[1]
        _dict['trembl_ec_number'] = i[2]
        _dict['trembl_go_component'] = i[3]
        _dict['trembl_go_process'] = i[4]
        _dict['trembl_go_function'] = i[5]
        _dict['trembl_interpro'] = i[6]
        _dict['trembl_pfam'] = i[7]
        records.append(tuple([_dict['trembl_name'], _dict['trembl_ec_number'], _dict['trembl_go_component'],
                              _dict['trembl_go_process'], _dict['trembl_go_function'], _dict['trembl_interpro'],
                              _dict['trembl_pfam'], _dict['trembl100_acc']]))
        record_counter += 1

    if record_counter > 0:
        print('Processing remaining tasks.')
        cursor_sql_db.execute("begin")
        cursor_sql_db.executemany(
            "UPDATE summary SET trembl_name = ?, trembl_ec_number = ?, trembl_go_component = ?,"
            "trembl_go_process = ?, trembl_go_function = ?, trembl_interpro = ?, trembl_pfam = ? where trembl100_acc=?",
            records)
        cursor_sql_db.execute("commit")
    cursor_sql_db.execute("CREATE INDEX summary_idx_trembl_dat ON summary (trembl_ec_number, trembl_go_component,"
                          "trembl_go_process, trembl_go_function, trembl_interpro, trembl_pfam)")

    time2 = time.time()
    time_consum = time2 - time1
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Finished parsing trembl dat_file in {:.2f} secs.'.format(time_consum))


def parse_pdb(pdb_output, sql_db):
    time1 = time.time()
    print('=' * 30)
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Starting to parse PDB output, please wait...')
    connect = sqlite3.connect(sql_db)
    cursor = connect.cursor()
    with open(pdb_output, 'r') as rd:
        text = [i.split('\t')[:2] for i in rd.readlines()]
    record_counter = 0
    records = []
    for i in text:
        if record_counter == 10000:
            cursor.execute("begin")
            cursor.executemany("UPDATE summary SET PDB100_acc = ? where acc_id=?", records)
            cursor.execute("commit")
            record_counter = 0
            records = []

        _dict = {}
        _dict['acc_id'] = i[0]
        _dict['pdb100_acc'] = i[1]
        records.append(tuple([_dict['pdb100_acc'], _dict['acc_id']]))
        record_counter += 1

    if record_counter > 0:
        print('Processing remaining tasks.')
        cursor.execute("begin")
        cursor.executemany("UPDATE summary SET PDB100_acc = ? where acc_id=?", records)
        cursor.execute("commit")
    cursor.execute("CREATE INDEX summary_idx_pdb ON summary (pdb100_acc)")

    time2 = time.time()
    time_consum = time2 - time1
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Finished parsing PDB output in {:.2f} secs.'.format(time_consum))


# should update the mmcif file into sql_db when crawling the data
def mmcif_sql(struct_path, sql_db):
    time1 = time.time()
    print('=' * 30)
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Starting to parse mmcif to sql, please wait...')
    connect = sqlite3.connect(sql_db)
    cursor = connect.cursor()
    record_counter = 0
    records = []
    struct_path_list = [os.path.join(struct_path, i) for i in os.listdir(struct_path)]
    for struct in struct_path_list:
        with open(struct, 'r') as rd:
            text = ''.join(rd.readlines())
        if record_counter == 1000:
            print('Saving 1000 mmcif files into sql, please wait...')
            cursor.execute("begin")
            cursor.executemany(f"UPDATE summary SET alphafold_mmcif_file = ? where acc_id=?", records)
            cursor.execute("commit")
            record_counter = 0
            records = []

        _dict = {}
        struct_file_name_pref = os.path.split(struct)[1].split('.')[0]
        _dict['acc_id'] = '_'.join(struct_file_name_pref.split('_')[:6])
        _dict['mmcif_file'] = text
        records.append(tuple([_dict['mmcif_file'], _dict['acc_id']]))
        record_counter += 1

    if record_counter > 0:
        print('Processing remaining tasks.')
        cursor.execute("begin")
        cursor.executemany(f"UPDATE summary SET alphafold_mmcif_file = ? where acc_id=?", records)
        cursor.execute("commit")
    print('Creating index...')
    cursor.execute("CREATE INDEX summary_idx_mmcif_file ON summary (alphafold_mmcif_file)")

    time2 = time.time()
    time_consum = time2 - time1
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Finished parse mmcif to sql in {:.2f} secs.'.format(time_consum))


def sele_mmcif(condition, sql_db, sele_str_dir_path):
    time1 = time.time()
    print('=' * 30)
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Starting to parse mmcif to sql, please wait...')
    connect = sqlite3.connect(sql_db)
    cursor = connect.cursor()
    sql = f"SELECT acc_id, cazy_family, alphafold_mmcif_file from summary where cazy_family like '{condition}%' " \
          f"and (sprot100_acc != '' or trembl100_acc != '')"
    cursor.execute(sql)
    ret_list = cursor.fetchall()
    for acc_id, cazy_family, af_file in ret_list:
        file_path = os.path.join(sele_str_dir_path, f'{acc_id}_' + cazy_family + '.cif')
        print('acc_id ', acc_id)
        if af_file:
            with open(file_path, 'w') as wt:
                wt.writelines(af_file)

    time2 = time.time()
    time_consum = time2 - time1
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Finished saving mmcif files in {:.2f} secs.'.format(time_consum))


def main():
    time1 = time.time()
    print('=' * 30)
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Starting to summary outputs, please wait...')

    cazy_output = r'D:\subject\active\PyMulstruct\data\tibet\cazy_overview.txt'
    merops_output = r'D:\subject\active\PyMulstruct\data\tibet\merops_output.tab'
    sql_db = r'D:\subject\active\PyMulstruct\data\tibet\summary.db'
    sprot_output = r'D:\subject\active\PyMulstruct\data\tibet\sprot_100.tab'
    sprot_db = r'G:\DB\uniprot_sprot.db'
    trembl_output = r'D:\subject\active\PyMulstruct\data\tibet\trembl_100.tab'
    trembl_db = r'G:\DB\newdb'
    pdb_output = r'D:\subject\active\PyMulstruct\data\tibet\pdb_100.tab'
    struct_path = r'D:\subject\active\PyMulstruct\data\tibet\str'
    sele_str_dir_path = r'D:\subject\active\PyMulstruct\data\tibet\target_str'

    create_sqlite(sql_db)
    parse_cazy(cazy_output, sql_db)
    parse_merops(merops_output, sql_db)
    parse_sprot_dmnd(sprot_output, sql_db)
    annot_sprot_sql(sprot_db, sql_db)
    parse_trembl_dmnd(trembl_output, sql_db)
    annot_trembl_sql(trembl_db, sql_db)
    parse_pdb(pdb_output, sql_db)
    # mmcif_sql(struct_path, sql_db)
    # sele_mmcif(condition='GH', sql_db=sql_db, sele_str_dir_path=sele_str_dir_path)

    time2 = time.time()
    time_consum = time2 - time1
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print('Totally finished summarizing in {:.2f} secs.'.format(time_consum))


main()
