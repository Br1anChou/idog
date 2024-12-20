#!/bin/bash -e

# @Usage      :bash fasterq-dump.sh <sra_dir> <out_dir> <sample_id> <ncpus>

###software vision###
#sratoolkit:3.0.2
#pigz:2.4

###parameter###
sra_dir=$1
out_dir=$2
sample_id=$3
ncpus=$4
##############

source activate RNAseq_E4

mkdir -p ${out_dir}/${sample_id}/reads
cd ${out_dir}/${sample_id}/reads

vdb-validate --version
vdb-validate ${sra_dir}/${sample_id}.sra

fasterq-dump \
-e ${ncpus} --split-3 \
--outdir ${out_dir}/${sample_id}/reads \
${sra_dir}/${sample_id}.sra

pigz -p ${ncpus} ./${sample_id}_1.fastq
pigz -p ${ncpus} ./${sample_id}_2.fastq
pigz -p ${ncpus} ./${sample_id}.fastq

#rm -rf ./${sample_id}_1.fastq
#rm -rf ./${sample_id}_2.fastq

conda deactivate
