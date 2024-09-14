#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :classifyGroups.py -h
# @Time       :2024/6/13 14:52
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :idog
# @Version    :python 3.10.6
# @Description:Write IDs to a separate sample list file based on grouping information(tsv format)
# @Usage      :python classifyGroups.py -h

import csv
import os

def process_tsv(input_file):
	# Create a dictionary to store the list of groups and corresponding run_ids
	group_dict = {}

	# read tsv file
	with open(input_file, newline='') as tsvfile:
		reader = csv.DictReader(tsvfile, delimiter='\t')
		for row in reader:
			run_id = row['run_id']
			group = row['group']
			if group not in group_dict:
				group_dict[group] = []
			group_dict[group].append(run_id)

	# Create a new file for each group and write in the corresponding run_id
	os.chdir(os.path.dirname(input_file))
	print("Current working directory:", os.getcwd())

	for group, run_ids in group_dict.items():
		with open(f"{group}.list", 'w') as outfile:
			for run_id in run_ids:
				outfile.write(f"{run_id}\n")

# example
process_tsv('/PATH/groups.txt')

