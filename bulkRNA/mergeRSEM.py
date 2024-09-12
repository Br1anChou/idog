#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :mergeRSEM
# @Time       :2024/6/6 16:42
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :idog
# @Version    :python 3.10.6
# @Description:
# @Usage      :

import os
import pandas as pd
import argparse

def merge_gene_expression(input_dir, sample_list_file, output_file, column_choice):
	# 映射用户输入的列名到实际列名
	column_map = {
		'TPM': 'TPM',
		'transcript_ids': 'transcript_id(s)',
		'length': 'length',
		'effective_length': 'effective_length',
		'expected_count': 'expected_count',
		'FPKM': 'FPKM'
	}

	if column_choice not in column_map:
		raise ValueError(f"Invalid column choice: {column_choice}. Must be one of {list(column_map.keys())}")

	actual_column = column_map[column_choice]

	# 读取样本列表文件
	if not os.path.exists(sample_list_file):
		raise FileNotFoundError(f"The file {sample_list_file} does not exist.")

	with open(sample_list_file, 'r') as f:
		sample_dirs = [line.strip() for line in f]

	# 检查是否有样本目录名
	if not sample_dirs:
		raise ValueError("The sample list is empty.")

	# 初始化一个空的DataFrame用于存储结果
	merged_df = pd.DataFrame()
	missing_samples = []  # 用于记录不存在的样本ID

	# 遍历每个样本目录
	for sample in sample_dirs:
		# 定义每个样本的文件路径
		file_path = os.path.join(input_dir, sample, "rsem", f"{sample}_rsem.genes.results")

		# 检查文件是否存在
		if not os.path.exists(file_path):
			print(f"Warning: The file {file_path} does not exist. Skipping sample {sample}.")
			missing_samples.append(sample)
			continue

		# 读取文件
		try:
			df = pd.read_csv(file_path, sep='\t', usecols=['gene_id', actual_column])
		except Exception as e:
			print(f"Error reading {file_path}: {e}. Skipping sample {sample}.")
			missing_samples.append(sample)
			continue

		# 重命名选择列为样本ID
		df = df.rename(columns={actual_column: sample})

		# 如果merged_df是空的，初始化它
		if merged_df.empty:
			merged_df = df
		else:
			# 合并数据
			merged_df = pd.merge(merged_df, df, on='gene_id', how='outer')

	# 保存结果到一个新的文件
	merged_df.to_csv(output_file, sep='\t', index=False)
	print(f"Merged gene expression matrix saved to {output_file}")

	# 输出不存在的样本ID
	if missing_samples:
		print("The following samples were not found or had errors and were skipped:")
		for sample in missing_samples:
			print(f"- {sample}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Merge gene expression matrices.')
	parser.add_argument('-i', '--input_dir', required=True, help='Directory containing sample directories.')
	parser.add_argument('-l', '--sample_list', required=True, help='File containing list of sample directories.')
	parser.add_argument('-o', '--output_file', required=True, help='Output file to save the merged matrix.')
	parser.add_argument('-c', '--column', required=True, choices=['TPM', 'transcript_ids', 'length', 'effective_length', 'expected_count', 'FPKM'], help='Column to extract from each file.')

	args = parser.parse_args()

	merge_gene_expression(args.input_dir, args.sample_list, args.output_file, args.column)
