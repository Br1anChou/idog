#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :mergeRSEM.py
# @Time       :2024/6/6 16:42
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :idog
# @Version    :python 3.10.6
# @Description:merge gene expression matrix for RSEM
# @Usage      :python mergeRSEM.py -h

import os
import pandas as pd
import argparse

def merge_gene_expression(input_dir, sample_list_file, output_file, column_choice):
	# map to actual column names
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

	if not os.path.exists(sample_list_file):
		raise FileNotFoundError(f"The file {sample_list_file} does not exist.")

	with open(sample_list_file, 'r') as f:
		sample_dirs = [line.strip() for line in f]

	# check sample name in the sample list
	if not sample_dirs:
		raise ValueError("The sample list is empty.")

	# initialise an empty DataFrame for storing results
	merged_df = pd.DataFrame()
	missing_samples = []  # Used to record non-existent sample IDs

	# Iterate over each sample catalogue
	for sample in sample_dirs:
		# file path for each sample
		file_path = os.path.join(input_dir, sample, "rsem", f"{sample}_rsem.genes.results")

		# check files
		if not os.path.exists(file_path):
			print(f"Warning: The file {file_path} does not exist. Skipping sample {sample}.")
			missing_samples.append(sample)
			continue

		# read file
		try:
			df = pd.read_csv(file_path, sep='\t', usecols=['gene_id', actual_column])
		except Exception as e:
			print(f"Error reading {file_path}: {e}. Skipping sample {sample}.")
			missing_samples.append(sample)
			continue

		# Rename the selected column as sample ID
		df = df.rename(columns={actual_column: sample})

		# initialise merged_df if it is empty
		if merged_df.empty:
			merged_df = df
		else:
			# merge data
			merged_df = pd.merge(merged_df, df, on='gene_id', how='outer')

	# save file
	merged_df.to_csv(output_file, sep='\t', index=False)
	print(f"Merged gene expression matrix saved to {output_file}")

	# Output non-existing sample IDs
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
