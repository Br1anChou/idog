#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :calcuStat.py
# @Time       :2024/6/11 10:00
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :idog
# @Version    :python 3.10.6
# @Description:Calculate maximum, minimum, median and mean values for each gene
# @Usage      :python calcuStat.py -h

import pandas as pd
import argparse

def calculate_gene_statistics(input_file, output_file):
	# read tsv file
	df = pd.read_csv(input_file, sep='\t')

	# keep gene_id column
	gene_ids = df['gene_id']

	# Calculate maximum, minimum, median and mean values for each gene
	expression_values = df.drop(columns=['gene_id'])
	max_values = expression_values.max(axis=1)
	min_values = expression_values.min(axis=1)
	median_values = expression_values.median(axis=1)
	mean_values = expression_values.mean(axis=1)
	std_values = expression_values.std(axis=1)

	# Calculation of CV values (standard deviation / mean)
	cv_values = std_values / mean_values

	# create dataframe
	stats_df = pd.DataFrame({
		'gene_id': gene_ids,
		'max': max_values,
		'min': min_values,
		'median': median_values,
		'mean': mean_values,
		'cv': cv_values
	})

	# output tsv file
	stats_df.to_csv(output_file, sep='\t', index=False)
	print(f"Gene statistics saved to {output_file}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Calculate statistics for gene expression matrix.')
	parser.add_argument('-i', '--input_file', required=True, help='Input file containing the merged gene expression matrix.')
	parser.add_argument('-o', '--output_file', required=True, help='Output file to save the gene statistics.')

	args = parser.parse_args()

	calculate_gene_statistics(args.input_file, args.output_file)
