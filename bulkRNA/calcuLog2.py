#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :calcuLog2.py
# @Time       :2024/6/11 09:41
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :idog
# @Version    :python 3.10.6
# @Description:Gene expression matrix normalisation
# @Usage      :python3 calcuLog2.py -h

import pandas as pd
import argparse
import numpy as np

def calculate_log2_expression(input_file, output_file):
	# read file
	df = pd.read_csv(input_file, sep='\t')

	# select gene_id
	gene_ids = df['gene_id']

	# calculate log2(x+1)
	expression_values = df.drop(columns=['gene_id'])
	log2_values = np.log2(expression_values + 1)

	# create DataFrame with gene_id and log2(x+1)
	log2_df = pd.DataFrame(log2_values, columns=expression_values.columns)
	log2_df.insert(0, 'gene_id', gene_ids)

	# save file
	log2_df.to_csv(output_file, sep='\t', index=False)
	print(f"log2(x+1) transformed matrix saved to {output_file}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Calculate log2(x+1) for gene expression matrix.')
	parser.add_argument('-i', '--input_file', required=True, help='Input file containing the merged gene expression matrix.')
	parser.add_argument('-o', '--output_file', required=True, help='Output file to save the log2(x+1) transformed matrix.')

	args = parser.parse_args()

	calculate_log2_expression(args.input_file, args.output_file)
