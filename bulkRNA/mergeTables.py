#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :mergeTables.py
# @Time       :2024/6/11 20:36
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :idog
# @Version    :python 3.10.6
# @Description:calculate statistics for different subgroups by expression matrices
# @Usage      :python mergeTables.py -h

import pandas as pd
import argparse

def parse_breed_tissue(breed_tissue):
	# extract breed and tissue from the string, where breed is the first part and tissue is the remaining part
	parts = breed_tissue.split('_')
	breed = parts[0]
	tissue = '_'.join(parts[1:])
	tissue = tissue.split('.')[0]  # remove file name extension
	return breed, tissue

def merge_statistics_and_expression(statistics_file, expression_file, breed_tissue, output_file):
	# read statistics and expression matrix files
	stats_df = pd.read_csv(statistics_file, sep='\t')
	expression_df = pd.read_csv(expression_file, sep='\t')

	# ensure that the gene_id columns of both DataFrames are aligned
	if not stats_df['gene_id'].equals(expression_df['gene_id']):
		raise ValueError("The gene_id columns in the two files do not match.")

	# create tpm_list column
	expression_values = expression_df.drop(columns=['gene_id'])
	tpm_list = expression_values.apply(lambda row: '[' + ','.join(row.astype(str)) + ']', axis=1)

	# add tpm_list columns to dataframe
	stats_df['tpm_list'] = tpm_list

	# If the breed or tissue parameter is supplied, the breed and tissue are parsed and added to the DataFrame.
	if breed_tissue:
		breed, tissue = parse_breed_tissue(breed_tissue)
		group = f"{breed}_{tissue}"
		stats_df['breed'] = breed
		stats_df['tissue'] = tissue
		stats_df['group'] = group

	# Save the results to a new tsv file
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
