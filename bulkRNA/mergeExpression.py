#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :mergeExpression.py
# @Time       :2024/6/20 09:52
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :RNAseq_pipeline
# @Version    :python 3.10.6
# @Description:merge gene expression matrix
# @Usage      :python mergeExpression.py -i “expression_matrices/sample*.tsv” -o combined_expression_matrix.tsv

import os
import pandas as pd
import argparse
import glob

def merge_expression_matrices(input_pattern, output_file):
	# 使用 glob 模块匹配输入文件路径
	input_files = glob.glob(input_pattern)

	# 检查是否找到输入文件
	if not input_files:
		print(f"No files matched the pattern: {input_pattern}")
		return

	# 初始化一个空的 DataFrame，用于存储合并的结果
	merged_df = None

	for file in input_files:
		# 读取文件
		df = pd.read_csv(file, sep='\t')

		# 如果 merged_df 为空，初始化它为当前的 DataFrame
		if merged_df is None:
			merged_df = df
		else:
			# 否则，将当前的 DataFrame 合并到 merged_df
			merged_df = pd.merge(merged_df, df, on='gene_id', how='outer')

	# 保存合并后的结果到输出文件，使用制表符分隔
	merged_df.to_csv(output_file, sep='\t', index=False)
	print(f"Combined expression matrix saved to {output_file}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Merge multiple gene expression matrices based on gene_id.')
	parser.add_argument('-i', '--input_pattern', required=True, help='Pattern to match input expression matrix files (e.g., "expression_matrices/*.tsv").')
	parser.add_argument('-o', '--output_file', required=True, help='Output file to save the combined expression matrix.')

	args = parser.parse_args()

	merge_expression_matrices(args.input_pattern, args.output_file)
