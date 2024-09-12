#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :calcuStat
# @Time       :2024/6/11 10:00
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :idog
# @Version    :python 3.10.6
# @Description:
# @Usage      :

import pandas as pd
import argparse

def calculate_gene_statistics(input_file, output_file):
	# 读取输入文件
	df = pd.read_csv(input_file, sep='\t')

	# 保留gene_id列
	gene_ids = df['gene_id']

	# 计算每个基因的最大值、最小值、中位数和平均值
	expression_values = df.drop(columns=['gene_id'])
	max_values = expression_values.max(axis=1)
	min_values = expression_values.min(axis=1)
	median_values = expression_values.median(axis=1)
	mean_values = expression_values.mean(axis=1)
	std_values = expression_values.std(axis=1)

	# 计算CV值（标准差 / 平均值）
	cv_values = std_values / mean_values

	# 创建包含统计值的DataFrame
	stats_df = pd.DataFrame({
		'gene_id': gene_ids,
		'max': max_values,
		'min': min_values,
		'median': median_values,
		'mean': mean_values,
		'cv': cv_values
	})

	# 保存结果到一个新的文件，使用制表符分隔
	stats_df.to_csv(output_file, sep='\t', index=False)
	print(f"Gene statistics saved to {output_file}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Calculate statistics for gene expression matrix.')
	parser.add_argument('-i', '--input_file', required=True, help='Input file containing the merged gene expression matrix.')
	parser.add_argument('-o', '--output_file', required=True, help='Output file to save the gene statistics.')

	args = parser.parse_args()

	calculate_gene_statistics(args.input_file, args.output_file)
