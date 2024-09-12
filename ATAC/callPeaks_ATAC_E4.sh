#!/usr/bin/env bash

# @File       :callPeaks_ATAC_E4
# @Time       :2024/9/1 15:33
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :ATAC_pipeline
# @Version    :
# @Description:
# @Usage      :

#in_file=2.pooled.tagAlign.gz
in_file=$1
out_dir=$2
sample_id=$3

source activate /p300s/zhaowm_group/tangbx/software/miniconda3/envs/ATAC_E4

macs2 callpeak \
-t ${in_file} \
-n ${sample_id} \
--outdir ${out_dir} \
-f BED -g 2353542698 -p 0.01 --shift -100 --extsize 200 \
--nomodel -B --SPMR --keep-dup all --call-summits

conda deactivate
