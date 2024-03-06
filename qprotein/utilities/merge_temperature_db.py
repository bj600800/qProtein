"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/10/10

# Description: 
# ------------------------------------------------------------------------------
"""
import csv
import os


def read_csv(file):
	all_temperature = []
	with open(file, "r") as f:
		csv_reader = csv.reader(f)
		for row in csv_reader:
			all_temperature.append(row)
	return all_temperature


def merge(first_db, second_db, merge_db):
	db1_name = [i[0] for i in first_db]
	for db2 in second_db:
		if db2[0].replace(" ", "_") not in db1_name:
			db2[0] = db2[0].replace(" ", "_")
			first_db.append(db2)
	with open(merge_db, "w", newline='') as file:
		writer = csv.writer(file, delimiter='\t')
		writer.writerows(first_db)
		
db_dir = r"D:\subject\active\1-qProtein\data\temperature_db"
first_db_path = os.path.join(db_dir, "organism_temperature.csv")
first_db = read_csv(first_db_path)
second_db_path = os.path.join(db_dir, "TEMPURA.csv")
second_db = read_csv(second_db_path)
merge_db = os.path.join(db_dir, "temperature.csv")
merge(first_db, second_db, merge_db)