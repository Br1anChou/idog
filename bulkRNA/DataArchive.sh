#!/usr/bin/env bash

# @File       :DataArchive.sh
# @Time       :2024/6/11 18:19
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :idog
# @Version    :1.0.0
# @Description:Extract and rename transcript and gene expression matrix files
# @Usage      :bash DataArchive.sh <sample.list> <file_main_directory> <file_output_directory>


# Check input parameters
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <sample.list> <file_main_directory> <file_output_directory>"
    exit 1
fi

# parameters
sample_list=$1
main_dir=$2
output_dir=$3

# Create the destination directory (if it doesn't exist)
mkdir -p "$output_dir"

# Read the ID of each sample in the sample.list file.
while IFS= read -r sample_id; do
    # Archive and rename kallisto_transcript expression matrix file
    kallisto_file="${main_dir}/${sample_id}/kallisto_transcript/abundance.tsv"
    if [ -f "$kallisto_file" ]; then
        cp "$kallisto_file" "${output_dir}/${sample_id}_transcripts_expression.tsv"
    else
        echo "File $kallisto_file does not exist."
    fi

    # Archive and rename rsem_gene expression matrix file
    rsem_file="${main_dir}/${sample_id}/rsem/${sample_id}_rsem.genes.results"
    if [ -f "$rsem_file" ]; then
        cp "$rsem_file" "${output_dir}/${sample_id}_genes_expression.tsv"
    else
        echo "File $rsem_file does not exist."
    fi
done < "$sample_list"
