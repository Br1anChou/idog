#!/bin/bash

###software vision###
#sratoolkit:3.0.2
#pigz:2.4

###parameter###
software_dir=/PATH/sratoolkit.3.0.2-centos_linux64/bin
out_dir=/PATH
sample_id=SRR13339272
ncpus=10
##############

mkdir -p ${out_dir}/reads
cd ${out_dir}/reads

${software_dir}/vdb-validate ./$sample_id

${software_dir}/fasterq-dump \
-e ${ncpus} --split-3 \
--outdir ${out_dir}/reads \
./${sample_id}

pigz -p ${ncpus} ./${sample_id}_1.fastq
pigz -p ${ncpus} ./${sample_id}_2.fastq

rm -rf ./${sample_id}_1.fastq
rm -rf ./${sample_id}_2.fastq
