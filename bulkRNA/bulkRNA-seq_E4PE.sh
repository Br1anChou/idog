#!/bin/bash -e

###software vision###
#fastp:0.23.2
#kallisto:0.46.0
#STAR:2.73a
#rsem:1.3.3


###parameter###
ref_dir=/p300s/zhaowm_group/tangbx/reference/UU_Cfam_GSD_1.0
#fastq_dir=/p300s/zhaowm_group/tangbx/idog/test
fastq_dir=$1
#out_dir=/p300s/zhaowm_group/tangbx/idog/test
out_dir=$2
#r1_file=SRR19614475_1.fastq.gz
#r2_file=SRR19614475_2.fastq.gz
#sample_id=SRR19614475
sample_id=$3
#ncpus=10
ncpus=$4
#ramGB=120
ramGB=$5
fwd_prob=0.5
#############

source activate RNAseq_E4

export PERL5LIB="/p300s/zhaowm_group/tangbx/software/miniconda3/envs/RNAseq_E4/lib/perl5/5.32"

mkdir -p ${out_dir}/${sample_id}/fastp

fastp -g -q 5 -u 50 -n 5 \
-i ${fastq_dir}/${sample_id}/reads/${sample_id}_1.fastq.gz \
-I ${fastq_dir}/${sample_id}/reads/${sample_id}_2.fastq.gz \
-o ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R1.fastq.gz \
-O ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R2.fastq.gz \
-j ${out_dir}/${sample_id}/fastp/${sample_id}_fastp.json \
-h ${out_dir}/${sample_id}/fastp/${sample_id}_fastp.html \
-R "${sample_id}_fastp_report"

mkdir -p ${out_dir}/${sample_id}/kallisto_gene 

kallisto quant \
-t ${ncpus} --fusion --plaintext \
-i ${ref_dir}/kallisto.index/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.111.gene.fa.idx \
-o ${out_dir}/${sample_id}/kallisto_gene \
${out_dir}/${sample_id}/fastp/${sample_id}.clean.R1.fastq.gz \
${out_dir}/${sample_id}/fastp/${sample_id}.clean.R2.fastq.gz

mkdir -p ${out_dir}/${sample_id}/kallisto_transcript

kallisto quant \
-t ${ncpus} --fusion --plaintext \
-i ${ref_dir}/kallisto.index/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.111.transcript.fa.idx \
-o ${out_dir}/${sample_id}/kallisto_transcript \
${out_dir}/${sample_id}/fastp/${sample_id}.clean.R1.fastq.gz \
${out_dir}/${sample_id}/fastp/${sample_id}.clean.R2.fastq.gz

##rsem

mkdir -p ${out_dir}/${sample_id}/star

STAR \
--genomeDir ${ref_dir}/star.index \
--readFilesIn ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R1.fastq.gz ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R2.fastq.gz \
--outFileNamePrefix ${out_dir}/${sample_id}/star/${sample_id}.gsd. \
--readFilesCommand zcat \
--runThreadN ${ncpus} \
--genomeLoad NoSharedMemory \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMheaderCommentFile COfile.txt \
--outSAMheaderHD @HD VN:1.4 SO:coordinate \
--outSAMunmapped Within \
--outFilterType BySJout \
--outSAMattributes NH HI AS NM MD \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--sjdbScore 1 \
--limitBAMsortRAM ${ramGB}000000000

STAR \
--runMode inputAlignmentsFromBAM \
--inputBAMfile ${out_dir}/${sample_id}/star/${sample_id}.gsd.Aligned.sortedByCoord.out.bam \
--outWigType bedGraph \
--outWigStrand Stranded \
--outWigReferencesPrefix chr \
--outFileNamePrefix ${out_dir}/${sample_id}/star/${sample_id}.gsd.

mkdir -p ${out_dir}/${sample_id}/rsem

rsem-calculate-expression --version

rsem-calculate-expression \
--bam \
--estimate-rspd \
--seed 12345 \
-p ${ncpus} \
--no-bam-output \
--forward-prob ${fwd_prob} \
--paired-end \
${out_dir}/${sample_id}/star/${sample_id}.gsd.Aligned.toTranscriptome.out.bam \
$ref_dir/rsem.index/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa \
${out_dir}/${sample_id}/rsem/${sample_id}_rsem

conda deactivate
