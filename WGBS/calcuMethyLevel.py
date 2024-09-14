#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :calcuMethyLevel
# @Time       :2024/9/4 15:43
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :RNAseq_pipeline
# @Version    :python 3.10.6
# @Description:calculate the methylation level within the regions
# @Usage      :python calcuMethyLevel.py -i <methylation_CpG/CHG/CHH_context.txt.gz> -o <methylation_level.txt> -b <regions.bed>

import os
from typing import List
# os.chdir('/Users/zhoubw/DataSpellProjects/idog/WGBS')

import subprocess
import pandas as pd
from intervaltree import IntervalTree
import argparse
import gzip


def calcu_methy(txt_file, bed_file, output_file):
	# Read bed files, extract gene information and build interval trees
	bed_columns = ["chrom", "start", "end", "gene_name", "score", "strand"]
	bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=bed_columns)

	# Initialise the interval tree that stores each chromosome
	chrom_intervals = {}

	# Iterating bed files to construct interval trees
	for _, row in bed_df.iterrows():
		chrom = str(row["chrom"])
		start = row["start"]
		end = row["end"]
		gene_name = row["gene_name"]
		gene_length = end - start

		if chrom not in chrom_intervals:
			chrom_intervals[chrom] = IntervalTree()

		# 插入基因区间到区间树中
		chrom_intervals[chrom][start:end] = {
			"gene_name": gene_name,
			"plus_count": 0,
			"minus_count": 0,
			"gene_length": gene_length
		}

	# Identify compressed or uncompressed files
	if txt_file.endswith(".gz"):
		file_open = gzip.open
	else:
		file_open = open

	# Quickly find the gene regions
	with file_open(txt_file, "rt") as f:
		line_number = 0
		for line in f:
			line_number += 1
			parts = line.strip().split("\t")

			# Make sure the file has 5 columns per line
			if len(parts) != 5:
				print(f"Skipping line {line_number}: expected 5 columns, got {len(parts)}. Content: {parts}")
				continue

			try:
				_, sign, chrom, position, _ = parts
				chrom = str(chrom)  # Make sure the chromosome number is string
				position = int(position)
			except ValueError as e:
				print(f"Error parsing line {line_number}: {e}. Content: {parts}")
				continue

			# Find the gene region in which the locus is located
			if chrom in chrom_intervals:
				overlapping_intervals = chrom_intervals[chrom][position]
				for interval in overlapping_intervals:
					if sign == "+":
						interval.data["plus_count"] += 1
					elif sign == "-":
						interval.data["minus_count"] += 1

	# write output
	with open(output_file, "w") as f:
		f.write(
			"gene_name\tregion_length\tmethylated_count\tunmethylated_count\tmethylated_percentage\tunmethylated_percentage\n")
		for chrom in chrom_intervals:
			for interval in chrom_intervals[chrom]:
				gene_name = interval.data["gene_name"]
				gene_length = interval.data["gene_length"]
				plus_count = interval.data["plus_count"]
				minus_count = interval.data["minus_count"]
				plus_percentage = (plus_count / gene_length) * 100
				minus_percentage = (minus_count / gene_length) * 100
				f.write(f"{gene_name}\t{gene_length}\t{plus_count}\t{minus_count}\t"
				        f"{plus_percentage:.2f}\t{minus_percentage:.2f}\n")

# Sort and remove duplicates
def run_shell_command(input_file, sorted_file):
	if input_file.endswith(".gz"):
		cmd = f"""
        pigz -dc {input_file} | \\
        tail -n +2 | \\
        sort -k3,3 -k4,4n | \\
        awk '!seen[$2,$3,$4,$5]++' > {sorted_file}
        """
	else:
		cmd = f"""
        tail -n +2 {input_file} | \\
        sort -k3,3 -k4,4n | \\
        awk '!seen[$2,$3,$4,$5]++' > {sorted_file}
        """

	subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')

def main():
	parser = argparse.ArgumentParser(description="Process methylation data against gene intervals in BED format.")
	parser.add_argument("-i", "--input", required=True, help="Input txt or txt.gz file with methylation data")
	parser.add_argument("-o", "--output", required=True, help="Output file name")
	parser.add_argument("-b", "--bed", required=True, help="BED file with gene intervals")
	args = parser.parse_args()

	# Create a temporary sorted file name
	sorted_txt_file = f"{args.input}.sorted.tmp"

	# First, run the shell command to process the input file
	run_shell_command(args.input, sorted_txt_file)

	# Now process the sorted file with the main logic
	calcu_methy(sorted_txt_file, args.bed, args.output)

	# Remove the intermediate sorted file
	os.remove(sorted_txt_file)


if __name__ == "__main__":
	main()
