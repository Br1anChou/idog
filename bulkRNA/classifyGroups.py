#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :classifyGroups
# @Time       :2024/6/13 14:52
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :idog
# @Version    :python 3.10.6
# @Description:
# @Usage      :

import csv
import os

def process_tsv(input_file):
	# 创建一个字典来存储group和对应的run_id列表
	group_dict = {}

	# 读取TSV文件
	with open(input_file, newline='') as tsvfile:
		reader = csv.DictReader(tsvfile, delimiter='\t')
		for row in reader:
			run_id = row['run_id']
			group = row['group']
			if group not in group_dict:
				group_dict[group] = []
			group_dict[group].append(run_id)

	# 为每个group创建一个新的文件，并写入对应的run_id
	os.chdir(os.path.dirname(input_file))
	print("Current working directory:", os.getcwd())

	for group, run_ids in group_dict.items():
		with open(f"{group}.list", 'w') as outfile:
			for run_id in run_ids:
				outfile.write(f"{run_id}\n")

# 示例：调用函数并传入TSV文件的路径
process_tsv('/Users/zhoubw/Documents/01 Work/02 Project/05 iDog/RNA/disease/standard_disease.txt')

