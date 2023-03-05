# -*- coding: utf-8 -*-
# @Time    : 2023/1/9 15:11
# @Author  : Zhixin Dou
"""
parse uniprot dat file to sqlite3
"""
import gzip
import sqlite3



def parse_uniprot_dat(dat_file, sql_db, table_name):
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print(f'Start create sql_db for {table_name}.')
    connect = sqlite3.connect(sql_db)
    cursor = connect.cursor()
    cursor.execute("DROP TABLE IF EXISTS uniprot_dat")
    cursor.execute(f"CREATE TABLE {table_name} "
                   "(accession TEXT, gene_name TEXT, ec_number TEXT, "
                   "go_component TEXT, go_process TEXT, go_function TEXT, "
                   "interpro TEXT, pfam TEXT)")

    with gzip.open(dat_file, 'rt') as uniprot_dat:
        accession = ''
        gene_name = ''
        ec_number = []
        go_component = []
        go_process = []
        go_function = []
        interpro = []
        pfam = []

        records = []
        record_counter = 0
        for line in uniprot_dat:
            _list = []
            if line.startswith("AC"):
                # GET accession number, the first one
                accession = line.split()[1].replace(';', '')

            elif line.startswith("DE   RecName") and "Full=" in line:
                gene_name = line.split("Full=")[1]
                gene_name = gene_name.split("{")[0].rstrip()
                gene_name = gene_name.lower()

            elif line.startswith("DE") and "EC=" in line:
                ec_code = line.strip().split()[1]
                ec_code = ec_code.replace(";", "")
                ec_code = ec_code.replace("EC=", "").lstrip()
                ec_number.append(ec_code)

            elif "DR   GO;" in line:
                if "; C:" in line:
                    # cellular component
                    go = '|'.join(line.rstrip().split('; ')[1:3]).rstrip('.')
                    go_component.append(go)
                elif "; F:" in line:
                    # molecular function
                    go = '|'.join(line.rstrip().split('; ')[1:3]).rstrip('.')
                    go_function.append(go)
                elif "; P:" in line:
                    # biological process
                    go = '|'.join(line.rstrip().split('; ')[1:3]).rstrip('.')
                    go_process.append(go)

            elif "DR   InterPro" in line:
                ipr = '|'.join(line.strip().split('; ')[1:3]).rstrip('.')
                interpro.append(ipr)

            elif "DR   Pfam" in line:
                pf = '|'.join(line.strip().split('; ')[1:3]).rstrip('.')
                pfam.append(pf)

            elif "//\n" in line:
                records.append(tuple([accession, gene_name, ';'.join(ec_number), ';'.join(go_component),
                                      ';'.join(go_process), ';'.join(go_function), ';'.join(interpro),
                                      ';'.join(pfam)]))
                record_counter += 1
                accession = ''
                gene_name = ''
                ec_number = []
                go_component = []
                go_process = []
                go_function = []
                interpro = []
                pfam = []

            if record_counter == 50000:
                cursor.execute("begin")
                cursor.executemany(f"INSERT INTO {table_name} VALUES(?, ?, ?, ?, ?, ?, ?, ?)", records)
                cursor.execute("commit")
                print('Insert 50000 items.')
                records = []
                record_counter = 0

        if record_counter > 0:
            print('Procession remaining records.')
            cursor.execute("begin")
            cursor.executemany(f"INSERT INTO {table_name} VALUES(?, ?, ?, ?, ?, ?, ?, ?)", records)
            cursor.execute("commit")
        print(f'Creating index for {table_name}.')
        cursor.execute(f"CREATE INDEX {table_name}_idx ON {table_name} (accession, gene_name, ec_number,"
                       "go_component, go_process, go_function, interpro, pfam)")

        print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
        print(f'Finished sql_db building for {table_name}.')


def main():
    # create uniprot_sprot dat db
    sprot_dat = r'G:\DB\uniprot_sprot.dat.gz'
    sprot_sql_db = r'G:\DB\uniprot_sprot.db'
    parse_uniprot_dat(sprot_dat, sprot_sql_db, 'uniprot_sprot')

main()
