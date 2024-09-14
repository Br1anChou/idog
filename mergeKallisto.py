#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :mergeKallisto.py
# @Time       :2023/5/24 17:54
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :RNAseq_pipeline
# @Version    :python 3.10.6
# @Description:Combine the results of kallisto into an expression matrix
# @Usage      :python mergeKallisto.py -i input_dir -o output_file -g gene_col -c expression_col -d


import os
import argparse
import pandas as pd

def merge_kallisto_data(sample_dir, output_file, gene_col=1, expression_col=5, del_zero=False):
    dfs = []
    # Iterate through the sample catalogue
    for d in os.scandir(sample_dir):
        # Processing of catalogue data only
        if d.is_dir():
            abundance_file_path = os.path.join(d.path, 'abundance.tsv')
            if os.path.exists(abundance_file_path):
                # read file
                try:
                    df = pd.read_csv(abundance_file_path, sep='\t', usecols=[gene_col, expression_col])
                    df.columns = ['target_id', d.name]
                    dfs.append(df)
                except Exception as e:
                    # print error
                    print(f"Error reading {abundance_file_path}: {e}")
            else:
                # print log
                print(f"File not found: {abundance_file_path}")

    # merge dataframe
    if dfs:
        merged_df = pd.concat(dfs, axis=1)
        # remove duplicate
        merged_df = merged_df.loc[:,~merged_df.columns.duplicated()]
        # remove rows with all zeros
        if del_zero:
            merged_df = merged_df.loc[merged_df.iloc[:, 1:].sum(axis=1) != 0]
        # output tsv
        merged_df.to_csv(output_file, sep='\t', index=False)
    else:
        print("No valid data found to merge.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge kallisto gene expression data into a matrix file.')
    parser.add_argument('-i', '--input_dir', type=str, required=True, help='Input directory containing sample directories.')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output file name.')
    parser.add_argument('-g', '--gene_col', type=int, default=1, help='Column number of target_id. Default is 1.')
    parser.add_argument('-c', '--expression_col', type=int, default=5, help='Column number of tpm. Default is 5.')
    parser.add_argument('-d', '--del_zero', action='store_true', help='Delete rows with zero expression in all samples.')
    args = parser.parse_args()
    merge_kallisto_data(args.input_dir, args.output_file, args.gene_col - 1, args.expression_col - 1, args.del_zero)
