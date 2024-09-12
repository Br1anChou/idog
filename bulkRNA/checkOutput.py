#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :checkOutput
# @Time       :2024/4/17 14:52
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :idog
# @Version    :python 3.10.6
# @Description:
# @Usage      :


import os
import re

def check_directories_exist(directory):
	subdirectories = ['fastp', 'kallisto_gene', 'kallisto_transcript', 'rsem', 'star']
	missing_directories = []
	for subdirectory in subdirectories:
		if not os.path.exists(os.path.join(directory, subdirectory)):
			missing_directories.append(subdirectory)
	return missing_directories

if __name__ == "__main__":
	fastp_directory = "/data_group/liuyanhu/liuyanhu2/zhoubw/proj/10.iDog/CFAMCANFSHP.YPi1956"  # 将此处替换为你的SRR目录的实际路径
	missing = check_directories_exist(fastp_directory)
	if missing:
		print("以下目录不存在：")
		for directory in missing:
			print(directory)
	else:
		print("所有目录都存在。")

import os
import json

def check_files_exist(directory):
	file1 = os.path.join(directory, "CFAMCANFSHP.YPi1956.clean.R1.fastq.gz")
	file2 = os.path.join(directory, "CFAMCANFSHP.YPi1956.clean.R2.fastq.gz")
	if not (os.path.exists(file1) and os.path.exists(file2)):
		return False
	return True

def extract_fastp_info(json_file):
	with open(json_file, 'r') as f:
		data = json.load(f)
		summary = data["summary"]
		after_filtering = summary["after_filtering"]
		total_reads = after_filtering["total_reads"]
		total_bases = after_filtering["total_bases"]
		q30_bases = after_filtering["q30_bases"]
		q30_rate = after_filtering["q30_rate"]
		gc_content = after_filtering["gc_content"]
		duplication_rate = data["duplication"]["rate"]
		return total_reads, total_bases, q30_bases, q30_rate, gc_content, duplication_rate

if __name__ == "__main__":
	fastp_directory = "/data_group/liuyanhu/liuyanhu2/zhoubw/proj/10.iDog/CFAMCANFSHP.YPi1956/fastp"  # 将此处替换为你的SRR目录的实际路径
	if check_files_exist(fastp_directory):
		json_file = os.path.join(fastp_directory, "CFAMCANFSHP.YPi1956_fastp.json")
		if os.path.exists(json_file):
			total_reads, total_bases, q30_bases, q30_rate, gc_content, duplication_rate = extract_fastp_info(json_file)
			print("total_reads:", total_reads)
			print("total_bases:", total_bases)
			print("q30_bases:", q30_bases)
			print("q30_rate:", q30_rate)
			print("gc_content:", gc_content)
			print("duplication_rate:", duplication_rate)
		else:
			print("找不到fastp的JSON文件。")
	else:
		print("CFAMCANFSHP.YPi1956.clean.R1.fastq.gz或CFAMCANFSHP.YPi1956.clean.R2.fastq.gz文件不存在。")


def check_bam_files_exist_and_nonempty(directory):
	bam_file_patterns = ['Aligned.sortedByCoord.out.bam', 'Aligned.toTranscriptome.out.bam']
	found_files = {pattern: False for pattern in bam_file_patterns}

	for file_name in os.listdir(directory):
		for pattern in bam_file_patterns:
			if file_name.endswith(pattern) and os.path.getsize(os.path.join(directory, file_name)) > 0:
				found_files[pattern] = True

	# Check if all required files are found and non-empty
	if all(found_files.values()):
		return True
	else:
		return False

def check_log_file(log_file):
	if not os.path.exists(log_file):
		return False

	with open(log_file, 'r') as f:
		log_content = f.read()
		# Updated regex to match dynamic date and time
		if re.search(r"\w{3} \d{2} \d{2}:\d{2}:\d{2} .* done", log_content):
			return True

	return False

if __name__ == "__main__":
	star_directory = "/data_group/liuyanhu/liuyanhu2/zhoubw/proj/10.iDog/CFAMCANFSHP.YPi1956/star"  # 将此处替换为你的SRR目录的实际路径
	if check_bam_files_exist_and_nonempty(star_directory):
		log_file = os.path.join(star_directory, "CFAMCANFSHP.YPi1956.gsd.Log.out")
		if check_log_file(log_file):
			print("STAR finished")
		else:
			print("STAR unfinished")
	else:
		print("bam file does not exist or file size is 0.")


import os

def check_abundance_file(directory):
	abundance_file = os.path.join(directory, "abundance.tsv")
	if os.path.exists(abundance_file):
		return abundance_file
	else:
		return None

def calculate_sum_of_columns(file_path, est_counts_index, tpm_index):
	with open(file_path, 'r') as file:
		next(file)  # skip header
		est_counts_sum = 0.0
		tpm_sum = 0.0
		for line in file:
			parts = line.strip().split('\t')
			est_counts_sum += float(parts[est_counts_index])
			tpm_sum += float(parts[tpm_index])
		return est_counts_sum, tpm_sum

if __name__ == "__main__":
	kallisto_gene_directory = "/data_group/liuyanhu/liuyanhu2/zhoubw/proj/10.iDog/CFAMCANFSHP.YPi1956/kallisto_gene"  # Set the actual path
	abundance_file = check_abundance_file(kallisto_gene_directory)

	if abundance_file:
		# Assuming est_counts is column 3 and tpm is column 4 (0-based index)
		est_counts_sum, tpm_sum = calculate_sum_of_columns(abundance_file, 3, 4)
		if est_counts_sum == 0.0 or tpm_sum == 0.0:
			print("计算结果错误")
		else:
			print("kallisto_gene finished!")
	else:
		print("文件不存在")

if __name__ == "__main__":
	kallisto_gene_directory = "/data_group/liuyanhu/liuyanhu2/zhoubw/proj/10.iDog/CFAMCANFSHP.YPi1956/kallisto_transcript"  # Set the actual path
	abundance_file = check_abundance_file(kallisto_gene_directory)

	if abundance_file:
		# Assuming est_counts is column 3 and tpm is column 4 (0-based index)
		est_counts_sum, tpm_sum = calculate_sum_of_columns(abundance_file, 3, 4)
		if est_counts_sum == 0.0 or tpm_sum == 0.0:
			print("计算结果错误")
		else:
			print("kallisto_transcript finished!")
	else:
		print("文件不存在")


import os

def find_files(directory, suffix):
	matching_files = []
	for filename in os.listdir(directory):
		if filename.endswith(suffix) and os.path.getsize(os.path.join(directory, filename)) > 0:
			matching_files.append(os.path.join(directory, filename))
	return matching_files

def check_sums(file_path, indices):
	with open(file_path, 'r') as file:
		header = next(file).strip().split('\t')
		sums = [0.0] * len(indices)
		for line in file:
			parts = line.strip().split('\t')
			for i, idx in enumerate(indices):
				sums[i] += float(parts[idx])
		return sums

if __name__ == "__main__":
	rsem_directory = "/data_group/liuyanhu/liuyanhu2/zhoubw/proj/10.iDog/CFAMCANFSHP.YPi1956/rsem"  # Set the actual path
	results = {}

	for suffix in ['_rsem.genes.results', '_rsem.isoforms.results']:
		files = find_files(rsem_directory, suffix)
		if not files:
			print(f"No valid files ending with {suffix}")
			continue

		for file_path in files:
			expected_count_idx = 4  # Assuming expected_count is the 5th column
			TPM_idx = 5           # Assuming TPM is the 6th column
			FPKM_idx = 6          # Assuming FPKM is the 7th column
			indices = [expected_count_idx, TPM_idx, FPKM_idx]
			sums = check_sums(file_path, indices)

			if all(s > 0 for s in sums):
				print(f"{os.path.basename(file_path)}: Data check passed!")
			else:
				print(f"{os.path.basename(file_path)}: Data may be all zeros.")
