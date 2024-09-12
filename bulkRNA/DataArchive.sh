#!/usr/bin/env bash

# @File       :DataArchive
# @Time       :2024/6/11 18:19
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :idog
# @Version    :
# @Description:
# @Usage      :

#!/bin/bash

# 检查是否提供了三个参数
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <sample.list> <file_main_directory> <file_output_directory>"
    exit 1
fi

# 读取传入的参数
sample_list=$1
main_dir=$2
output_dir=$3

# 创建目标目录（如果不存在的话）
mkdir -p "$output_dir"

# 读取 sample.list 文件中的每个样本 ID
while IFS= read -r sample_id; do
    # 归档并重命名 kallisto_transcript 文件
    kallisto_file="${main_dir}/${sample_id}/kallisto_transcript/abundance.tsv"
    if [ -f "$kallisto_file" ]; then
        cp "$kallisto_file" "${output_dir}/${sample_id}_transcripts_expression.tsv"
    else
        echo "File $kallisto_file does not exist."
    fi

    # 归档并重命名 rsem 文件
    rsem_file="${main_dir}/${sample_id}/rsem/${sample_id}_rsem.genes.results"
    if [ -f "$rsem_file" ]; then
        cp "$rsem_file" "${output_dir}/${sample_id}_genes_expression.tsv"
    else
        echo "File $rsem_file does not exist."
    fi
done < "$sample_list"
