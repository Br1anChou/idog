#!/usr/bin/env bash

# @File       :ATAC_E4.sh
# @Time       :2024/8/23 10:44
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :ATAC_pipeline
# @Version    :1.0.0
# @Description:main script for ATAC_pipeline
# @Usage      :bash ATAC_E4.sh <fastq_dir> <out_dir> <sample_id> <ncpus> <ramGB>

###software vision###
#fastp:0.23.2
#bowte2:2.5.3
#samtools:1.15.1
#picard:2.22.8
#bedtools:2.31.1
#macs2:2.2.4

###parameter###
ref_dir=/path/reference/UU_Cfam_GSD_1.0
tmp_dir=/p300s/zhaowm_group/tangbx/idog/tmp
#fastq_dir=/gpfs/zhaowm_group/wangyibo/epigenome/ATAC/PRJEB20256
fastq_dir=$1
#out_dir=/p300s/zhaowm_group/tangbx/idog/atac_result
out_dir=$2
#r1_file=SRR8090327.fastq.gz
#sample_id=ERR1912239
sample_id=$3
#ncpus=10
ncpus=$4
#ramGB=60
ramGB=$5
#############

###private environment###
script=/p300s/zhaowm_group/tangbx/idog/test/atac
software=/software/biosoft/software/picard2.22.8/picard/build/libs
source activate /p300s/zhaowm_group/tangbx/software/miniconda3/envs/RNAseq_E4
export PERL5LIB="/p300s/zhaowm_group/tangbx/software/miniconda3/envs/RNAseq_E4/lib/perl5/5.32"
export JAVA_HOME=/software/biosoft/software/jdk1.8/jdk1.8.0_45
export PATH=$JAVA_HOME/bin:$JAVA_HOME/jre/bin:$PATH
export CLASSPATH=$JAVA_HOME/lib:$JAVA_HOME/jre/lib:$CLASSPATH
export JRE_HOME=$JAVA_HOME/jre
############

###fastp###
mkdir -p ${out_dir}/${sample_id}/fastp
# 检查是否存在以_1.fastq.gz或_2.fastq.gz结尾的文件
PE_raw=$(find "${fastq_dir}/${sample_id}/reads/" -type f \( -name "${sample_id}_1.fastq.gz" -o -name "${sample_id}_2.fastq.gz" \))
# 检查是否存在以.fastq.gz结尾的文件
SE_raw=$(find "${fastq_dir}/${sample_id}/reads/" -type f -name "${sample_id}.fastq.gz")

###判断条件并执行相应的代码
if [[ -n "$PE_raw" && -z "$SE_raw" ]]; then
    # PE代码的内容
    echo "...run fastp PE pipeline..."
    paired_end=0
    fastp -g -q 5 -u 50 -n 5 \
    -i ${fastq_dir}/${sample_id}/reads/${sample_id}_1.fastq.gz \
    -I ${fastq_dir}/${sample_id}/reads/${sample_id}_2.fastq.gz \
    -o ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R1.fastq.gz \
    -O ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R2.fastq.gz \
    -j ${out_dir}/${sample_id}/fastp/${sample_id}_fastp.json \
    -h ${out_dir}/${sample_id}/fastp/${sample_id}_fastp.html \
    -R "${sample_id}_fastp_report"

    cd ${out_dir}/${sample_id}/fastp
    json_file="${sample_id}_fastp.json"

    read1_length=$(jq -r '.summary.after_filtering.read1_mean_length' "$json_file")
    read2_length=$(jq -r '.summary.after_filtering.read2_mean_length' "$json_file")

    echo "read1_mean_length: $read1_length"
    echo "read2_mean_length: $read2_length"

elif [[ -z "$PE_raw" && -n "$SE_raw" ]]; then
    # SE代码的内容
    echo "...run fastp SE pipeline..."
    paired_end=1
    fastp -g -q 5 -u 50 -n 5 \
    -i ${fastq_dir}/${sample_id}/reads/${sample_id}.fastq.gz \
    -o ${out_dir}/${sample_id}/fastp/${sample_id}.clean.fastq.gz \
    -j ${out_dir}/${sample_id}/fastp/${sample_id}_fastp.json \
    -h ${out_dir}/${sample_id}/fastp/${sample_id}_fastp.html \
    -R "${sample_id}_fastp_report"

    cd ${out_dir}/${sample_id}/fastp
    json_file="${sample_id}_fastp.json"

    read1_length=$(jq -r '.summary.after_filtering.read1_mean_length' "$json_file")
    echo "read1_mean_length: $read1_length"

else
    echo "ERROR: Raw sequencing file dose not exist or have a wrong name!"
    exit 1
fi


###bowtie2###
mkdir -p ${out_dir}/${sample_id}/bowtie2
cd ${out_dir}/${sample_id}/bowtie2

if [ "$paired_end" -eq 0 ]; then
  echo "...run bowtie2 PE pipeline..."
  bowtie2 --mm --threads ${ncpus} -X2000 -q \
  -x $ref_dir/bowtie2.index/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa \
  -1 ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R1.fastq.gz \
  -2 ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R2.fastq.gz | \
  samtools view -@ ${ncpus} -1 -S -b > ${sample_id}.raw.bam

  samtools view -@ ${ncpus} -F 1804 -q 30 -f 2 -u ${sample_id}.raw.bam | \
  samtools sort -@ ${ncpus} -m 2G -n -O bam -T ${tmp_dir} -o ${sample_id}.tmp.bam
  samtools fixmate -@ ${ncpus} -r -O bam ${sample_id}.tmp.bam ${sample_id}.fixmate.bam

  samtools view -@ ${ncpus} -F 1804 -f 2 -u ${sample_id}.fixmate.bam |
  samtools sort -@ ${ncpus} -m 2G -O bam -T ${tmp_dir} -o ${sample_id}.filt.bam

elif [ "$paired_end" -eq 1 ]; then
  echo "...run bowtie2 SE pipeline..."
  bowtie2 --mm --threads ${ncpus} \
  -x $ref_dir/bowtie2.index/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa \
  -U ${out_dir}/${sample_id}/fastp/${sample_id}.clean.fastq.gz | \
  samtools view -@ ${ncpus} -1 -S -b > ${sample_id}.raw.bam

  samtools view -@ ${ncpus} -F 1804 -q 30 -u ${sample_id}.raw.bam | \
  samtools sort -@ ${ncpus} -m 2G -n -O bam -T ${tmp_dir} -o ${sample_id}.filt.bam

else
    echo "ERROR: ${sample_id}.clean.R1|R2.fastq.gz does not exist!"
    exit 1
fi

samtools index -@ ${ncpus} -b ${sample_id}.filt.bam

java -Xmx4G -XX:ParallelGCThreads=1 \
-jar ${software}/picard.jar MarkDuplicates \
INPUT=${sample_id}.filt.bam \
OUTPUT=${sample_id}.dupmark.bam \
METRICS_FILE=${sample_id}.dup.qc \
VALIDATION_STRINGENCY=LENIENT \
USE_JDK_DEFLATER=TRUE \
USE_JDK_INFLATER=TRUE \
ASSUME_SORTED=TRUE \
REMOVE_DUPLICATES=FALSE

if [ "$paired_end" -eq 0 ]; then
  echo "...run remove duplication PE pipeline..."
  samtools view -@ ${ncpus} -F 1804 -f 2 -b ${sample_id}.dupmark.bam > ${sample_id}.nodup.bam
elif [ "$paired_end" -eq 1 ]; then
  echo "...run remove duplication SE pipeline..."
  samtools view -@ ${ncpus} -F 1804 -b ${sample_id}.dupmark.bam > ${sample_id}.nodup.bam
else
  echo "ERROR: ${sample_id}.dupmark.bam does not exist!"
  exit 1
fi

samtools index -@ ${ncpus} -b ${sample_id}.nodup.bam
samtools stats -@ ${ncpus} ${sample_id}.nodup.bam > ${sample_id}.nodup.bam.stats
samtools flagstat -@ ${ncpus} -O tsv ${sample_id}.nodup.bam > ${sample_id}.nodup.bam.flagstat

#rm -f ${sample_id}.sorted.bam

###TAG-Align###
if [ "$paired_end" -eq 0 ]; then
  echo "...run bam2ta PE pipeline..."
  bedtools bamtobed -bedpe -mate1 -i ${sample_id}.nodup.bam | \
  awk 'BEGIN{OFS="\t"} {printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n", $1,$2,$3,$9,$4,$5,$6,$10}' | \
  gzip -nc > ${sample_id}.nodup.tagAlign.gz

elif [ "$paired_end" -eq 1 ]; then
  echo "...run bam2ta SE pipeline..."
  bedtools bamtobed -i ${sample_id}.nodup.bam | \
  awk 'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print $0}}\' | \
  gzip -nc > ${sample_id}.nodup.tagAlign.gz

else
    echo "ERROR: ${sample_id}.nodup.bam does not exist!"
    exit 1
fi

zcat -f ${sample_id}.nodup.tagAlign.gz | awk 'BEGIN {OFS = "\t"} {
    if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5}
    if ($2 >= $3) {if ($6 == "+") {$2 = $3 - 1} else {$3 = $2 + 1}}
    print $0
}' | gzip -nc > ${sample_id}.tn5.tagAlign.gz

