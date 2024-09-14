#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :checkOutput_E4
# @Time       :2024/4/17 14:52
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :idog
# @Version    :python 3.10.6
# @Description:The input parameter path cannot include the sample ID
# @Usage      :python checkOutput_E4.py -p /results_path -i SRRXXXX


import os
import re
import json
import argparse

def file_exists_and_not_empty(file_path):
	"""Check if a file exists and is not empty."""
	return os.path.exists(file_path) and os.path.getsize(file_path) > 0

def sum_columns(file_path, indices):
	"""Sum specified columns in a tab-separated file."""
	sums = [0.0] * len(indices)
	with open(file_path, 'r') as file:
		next(file)  # skip header
		for line in file:
			parts = line.strip().split('\t')
			for i, idx in enumerate(indices):
				sums[i] += float(parts[idx])
	return sums

def main(base_path, sample_id):
	sample_id = args.id
	base_path = os.path.join(base_path, sample_id)

	print(f'Sample ID: {sample_id}')
	print(f'Base path: {base_path}')

	# Step 1: Check directories
	required_dirs = ["fastp", "kallisto_gene", "kallisto_transcript", "rsem", "star"]
	for d in required_dirs:
		dir_path = os.path.join(base_path, d)
		if not os.path.exists(dir_path):
			print(f"ERROR: {dir_path} missing!")
			return
		else:
			print(f"Directory exist: {dir_path}")


	# Step 2: Check fastp files and read JSON
	fastp_path = os.path.join(base_path, "fastp")
	json_file = os.path.join(fastp_path, sample_id+"_fastp.json")
	if file_exists_and_not_empty(json_file):
		with open(json_file, 'r') as jf:
			data = json.load(jf)
			after_filtering = data['summary']['after_filtering']
			# print(json.dumps(after_filtering, indent=4))
	else:
		print(f"ERROR: JSON file {json_file} missing or empty!")

	# Step 3: Check STAR outputs
	star_path = os.path.join(base_path, "star")
	log_file = os.path.join(star_path, sample_id+".gsd.Log.out")
	# check bam files
	# ref='gsd'
	bam_files = [sample_id+".gsd.Aligned.sortedByCoord.out.bam", sample_id+".gsd.Aligned.toTranscriptome.out.bam"]
	if all(file_exists_and_not_empty(os.path.join(star_path, f)) for f in bam_files):
		if file_exists_and_not_empty(log_file):
			with open(log_file) as lf:
				log_content = lf.read()
				# get finished time
				match = re.search(r"(\w{3} \d{2} \d{2}:\d{2}:\d{2}).*done", log_content)
				if match:
					finished_time = match.group(1)
					print(f"STAR finished: Completion date {finished_time}")
				else:
					print(f"ERROR: STAR log file {log_file} does not indicate completion.")
		else:
			print(f"ERROR: STAR log file {log_file} is missing or empty.")
	else:
		print(f"ERROR: STAR bam files {','.join(bam_files)} are missing or empty.")

	# Step 4: Check kallisto results
	kallisto_gene = os.path.join("kallisto_gene", "abundance.tsv")
	kallisto_gene_path = os.path.join(base_path, kallisto_gene)
	if file_exists_and_not_empty(kallisto_gene_path):
		if all(s > 0 for s in sum_columns(kallisto_gene_path, [3, 4])):
			print(f"kallisto finished: Data of {kallisto_gene} check passed!")
		else:
			print(f"ERROR: Data of {kallisto_gene} may be all zeros.")
	else:
		print(f"ERROR: {kallisto_gene_path} missing!")

	kallisto_transcript = os.path.join("kallisto_transcript", "abundance.tsv")
	kallisto_transcript_path = os.path.join(base_path, kallisto_transcript)
	if file_exists_and_not_empty(kallisto_transcript_path):
		if all(s > 0 for s in sum_columns(kallisto_transcript_path, [3, 4])):
			print(f"kallisto finished: Data of {kallisto_transcript} check passed!")
		else:
			print(f"ERROR: Data of {kallisto_transcript} may be all zeros.")
	else:
		print(f"ERROR: {kallisto_transcript_path} missing!")

	# Step 5: Check RSEM results
	rsem_path = os.path.join(base_path, "rsem")
	for suffix in ['_rsem.genes.results', '_rsem.isoforms.results']:
		files = [os.path.join(rsem_path, f) for f in os.listdir(rsem_path) if f.endswith(suffix) and file_exists_and_not_empty(os.path.join(rsem_path, f))]
		if not files:
			print(f"ERROR: {suffix} missing!")
			continue
		for file_path in files:
			if all(s > 0 for s in sum_columns(file_path, [4, 5, 6])):
				print(f"rsem finished: Data of {os.path.basename(file_path)} check passed!")
			else:
				print(f"ERROR: Data of {os.path.basename(file_path)} may be all zeros.")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Check analysis outputs for bulkRNA-seq_E4.sh pipeline.")
	parser.add_argument('-p', '--path', required=True, type=str, help='Base path where analysis directories are located')
	parser.add_argument('-i', '--id', required=True, type=str, help='Sample ID to check')

	args = parser.parse_args()

	main(args.path, args.id)
