#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :mergeTables.py
# @Time       :2024/6/11 20:36
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :idog
# @Version    :python 3.10.6
# @Description:
# @Usage      :python mergeTables.py -h

import pandas as pd
import argparse

def parse_breed_tissue(breed_tissue):
	# 从字符串中提取breed和tissue，breed是第一个部分，tissue是剩余部分
	parts = breed_tissue.split('_')
	breed = parts[0]
	tissue = '_'.join(parts[1:])
	tissue = tissue.split('.')[0]  # 去掉文件名的扩展名
	return breed, tissue

def merge_statistics_and_expression(statistics_file, expression_file, breed_tissue, output_file):
	# 读取统计文件和表达矩阵文件
	stats_df = pd.read_csv(statistics_file, sep='\t')
	expression_df = pd.read_csv(expression_file, sep='\t')

	# 确保两个DataFrame的gene_id列是对齐的
	if not stats_df['gene_id'].equals(expression_df['gene_id']):
		raise ValueError("The gene_id columns in the two files do not match.")

	# 创建TPM列表列
	expression_values = expression_df.drop(columns=['gene_id'])
	tpm_list = expression_values.apply(lambda row: '[' + ','.join(row.astype(str)) + ']', axis=1)

	# 添加TPM列表列到统计DataFrame中
	stats_df['tpm_list'] = tpm_list

	# 如果提供了breed_tissue参数，解析breed和tissue并添加到DataFrame中
	if breed_tissue:
		breed, tissue = parse_breed_tissue(breed_tissue)
		group = f"{breed}_{tissue}"
		stats_df['breed'] = breed
		stats_df['tissue'] = tissue
		stats_df['group'] = group

	# 保存结果到一个新的文件，使用制表符分隔
	stats_df.to_csv(output_file, sep='\t', index=False)
	print(f"Combined statistics and expression matrix saved to {output_file}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Merge gene statistics with expression matrix.')
	parser.add_argument('-s', '--statistics_file', required=True, help='Input file containing the gene statistics.')
	parser.add_argument('-e', '--expression_file', required=True, help='Input file containing the merged gene expression matrix.')
	parser.add_argument('-b', '--breed_tissue', required=False, help='String in the format breed_tissue.list to extract breed and tissue information.')
	parser.add_argument('-o', '--output_file', required=True, help='Output file to save the combined result.')

	args = parser.parse_args()

	merge_statistics_and_expression(args.statistics_file, args.expression_file, args.breed_tissue, args.output_file)
