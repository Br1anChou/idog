#!/usr/bin/env bash

# @File       :callPeaks_ATAC_E4.sh
# @Time       :2024/9/1 15:33
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :ATAC_pipeline
# @Version    :1.0.0
# @Description:ATAC peaks calling
# @Usage      :bash callPeaks_ATAC_E4.sh <in_file> <out_dir> <sample_id>

#in_file=sample.pooled.tagAlign.gz
in_file=$1
out_dir=$2
sample_id=$3

source activate /p300s/zhaowm_group/tangbx/software/miniconda3/envs/ATAC_E4

#dog genome size: 2353542698
macs2 callpeak \
-t ${in_file} \
-n ${sample_id} \
--outdir ${out_dir} \
-f BED -g 2353542698 -p 0.01 --shift -100 --extsize 200 \
--nomodel -B --SPMR --keep-dup all --call-summits

conda deactivate
