#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :calcuZscore.py
# @Time       :2024/6/7 12:33
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :idog
# @Version    :python 3.10.6
# @Description:calculate z-score
# @Usage      :python calcuZscore.py -h

import pandas as pd
import argparse
from scipy.stats import zscore

def calculate_zscore(input_file, output_file):
	# read file
	df = pd.read_csv(input_file, sep='\t')

	# select gene_id
	gene_ids = df['gene_id']

	# calculate z-score
	expression_values = df.drop(columns=['gene_id'])
	z_scores = expression_values.apply(zscore, axis=1)

	# create DataFrame with gene_id and z-scores
	zscore_df = pd.DataFrame(z_scores, columns=expression_values.columns)
	zscore_df.insert(0, 'gene_id', gene_ids)

	# save file
	zscore_df.to_csv(output_file, sep='\t', index=False)
	print(f"Z-scores matrix saved to {output_file}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Calculate z-scores for gene expression matrix.')
	parser.add_argument('-i', '--input_file', required=True, help='Input file containing the merged gene expression matrix.')
	parser.add_argument('-o', '--output_file', required=True, help='Output file to save the z-scores matrix.')

	args = parser.parse_args()

	calculate_zscore(args.input_file, args.output_file)
