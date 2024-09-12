#!/usr/bin/env bash -e

# @File       :WGBS_E3.sh
# @Time       :2024/7/11 16:40
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :WGBS_pipeline
# @Version    :1.0.0
# @Description:main script for BS-seq pipeline
# @Usage      :bash WGBS_E3.sh <fastq_dir> <out_dir> <sample_id> <ncpus> <ramGB>


###software vision###
#fastp:0.23.2
#bismark 0.23.1
#bowtie2:2.5.4

###parameter###
ref_dir=/p300s/zhaowm_group/tangbx/reference/UU_Cfam_GSD_1.0
tmp_dir=/p300s/zhaowm_group/tangbx/idog/tmp
software=/p300s/zhaowm_group/tangbx/idog/software/cgmaptools-0.1.3
#fastq_dir=/p300s/zhaowm_group/tangbx/idog/epi_results
fastq_dir=$1
#out_dir=/p300s/zhaowm_group/tangbx/idog/epi_results
out_dir=$2
#r1_file=SRR8090327.fastq.gz
#sample_id=SRR4064664
sample_id=$3
#ncpus=10
ncpus=$4
#ramGB=60
ramGB=$5
#############

###private environment###
source activate WGBS_E4
export PERL5LIB=/p300s/zhaowm_group/tangbx/software/miniconda3/envs/WGBS_E4/lib/perl5
# export PERL5LIB="/p300s/zhaowm_group/tangbx/software/miniconda3/envs/RNAseq_E4/lib/perl5/5.32"
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
    echo "...run PE pipeline..."
    paired_end=0
    fastp -g -q 5 -u 50 -n 5 \
    -i ${fastq_dir}/${sample_id}/reads/${sample_id}_1.fastq.gz \
    -I ${fastq_dir}/${sample_id}/reads/${sample_id}_2.fastq.gz \
    -o ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R1.fastq.gz \
    -O ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R2.fastq.gz \
    -j ${out_dir}/${sample_id}/fastp/${sample_id}_fastp.json \
    -h ${out_dir}/${sample_id}/fastp/${sample_id}_fastp.html \
    -R "${sample_id}_fastp_report"

elif [[ -z "$PE_raw" && -n "$SE_raw" ]]; then
    # SE代码的内容
    echo "...run SE pipeline..."
    paired_end=1
    fastp -g -q 5 -u 50 -n 5 \
    -i ${fastq_dir}/${sample_id}/reads/${sample_id}.fastq.gz \
    -o ${out_dir}/${sample_id}/fastp/${sample_id}.clean.fastq.gz \
    -j ${out_dir}/${sample_id}/fastp/${sample_id}_fastp.json \
    -h ${out_dir}/${sample_id}/fastp/${sample_id}_fastp.html \
    -R "${sample_id}_fastp_report"

else
    echo "ERROR: Raw sequencing file dose not exist or have a wrong name!"
    exit 1
fi

###run bismark
rm -rf ${out_dir}/${sample_id}/bismark
mkdir -p ${out_dir}/${sample_id}/bismark

# 检查是否存在以clean.R1.fastq.gz或clean.R2.fastq.gz结尾的文件
PE_clean=$(find "${out_dir}/${sample_id}/fastp/" -type f \( -name "${sample_id}.clean.R1.fastq.gz" -o -name "${sample_id}.clean.R2.fastq.gz" \))
# 检查是否存在以clean.fastq.gz结尾的文件
SE_clean=$(find "${out_dir}/${sample_id}/fastp/" -type f -name "${sample_id}.clean.fastq.gz")

if [[ -n "$PE_clean" && -z "$SE_clean" && "$paired_end" -eq 0 ]]; then
    # PE代码的内容
    echo "Input files: $PE_clean"
    bismark \
    --bowtie2 \
    --parallel $ncpus \
    --bam \
    --nucleotide_coverage \
    --genome_folder $ref_dir/bismark.index \
    --temp_dir ${tmp_dir}/${sample_id} \
    --path_to_bowtie2 /p300s/zhaowm_group/tangbx/software/miniconda3/envs/WGBS_E4/bin \
    -o ${out_dir}/${sample_id}/bismark \
    -1 ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R1.fastq.gz \
    -2 ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R2.fastq.gz

    deduplicate_bismark \
    --paired \
    --bam \
    --samtools_path /p300s/zhaowm_group/tangbx/software/miniconda3/envs/RNAseq_E4/bin/ \
    --output_dir ${out_dir}/${sample_id}/bismark \
    ${out_dir}/${sample_id}/bismark/${sample_id}.clean.R1_bismark_bt2_pe.bam

    bismark_methylation_extractor \
    --paired-end \
    --no_overlap \
    -o ${out_dir}/${sample_id}/bismark \
    --comprehensive \
    --parallel ${ncpus} \
    --buffer_size ${ramGB}G \
    --report \
    --bedGraph \
    --cytosine_report \
    --CX_context \
    --samtools_path /p300s/zhaowm_group/tangbx/software/miniconda3/envs/RNAseq_E4/bin/ \
    --genome_folder $ref_dir/bismark.index \
    ${out_dir}/${sample_id}/bismark/${sample_id}.clean.R1_bismark_bt2_pe.deduplicated.bam

elif [[ -z "$PE_clean" && -n "$SE_clean" && "$paired_end" -eq 1 ]]; then
    # SE代码的内容
    echo "Input files: $SE_clean"
    bismark \
    --bowtie2 \
    --parallel $ncpus \
    --bam \
    --nucleotide_coverage \
    --genome_folder $ref_dir/bismark.index \
    --temp_dir ${tmp_dir}/${sample_id} \
    --path_to_bowtie2 /p300s/zhaowm_group/tangbx/software/miniconda3/envs/WGBS_E4/bin \
    -o ${out_dir}/${sample_id}/bismark \
    ${out_dir}/${sample_id}/fastp/${sample_id}.clean.fastq.gz

    deduplicate_bismark \
    -s \
    --bam \
    --samtools_path /p300s/zhaowm_group/tangbx/software/miniconda3/envs/RNAseq_E4/bin/ \
    --output_dir ${out_dir}/${sample_id}/bismark \
    ${out_dir}/${sample_id}/bismark/${sample_id}.clean_bismark_bt2.bam

    bismark_methylation_extractor \
    -s \
    --no_overlap \
    -o ${out_dir}/${sample_id}/bismark \
    --comprehensive \
    --parallel ${ncpus} \
    --buffer_size ${ramGB}G \
    --report \
    --bedGraph \
    --cytosine_report \
    --CX_context \
    --samtools_path /p300s/zhaowm_group/tangbx/software/miniconda3/envs/RNAseq_E4/bin/ \
    --genome_folder $ref_dir/bismark.index \
    ${out_dir}/${sample_id}/bismark/${sample_id}.clean_bismark_bt2.deduplicated.bam

else
    echo "ERROR: Clean sequencing files does not exist or has a wrong name!"
    exit 1
fi

cd ${out_dir}/${sample_id}/bismark

bismark2report \
-o ${sample_id}_bismark_report \
--dir ${out_dir}/${sample_id}/bismark

if [ "$paired_end" -eq 0 ]; then
  bam_file="${out_dir}/${sample_id}/bismark/${sample_id}.clean.R1_bismark_bt2_pe.deduplicated.bam"
  CX_file="${sample_id}.clean.R1_bismark_bt2_pe.deduplicated.CX_report.txt"
  echo "Input bam files: $bam_file"
  echo "Input CX files: $CX_file"

  bam2nuc \
  --genome_folder $ref_dir/bismark.index \
  --dir ${out_dir}/${sample_id}/bismark \
  --samtools_path /p300s/zhaowm_group/tangbx/software/miniconda3/envs/RNAseq_E4/bin/ \
  ${out_dir}/${sample_id}/bismark/${sample_id}.clean.R1_bismark_bt2_pe.deduplicated.bam

  $software/cgmaptools convert bismark2cgmap \
  -i ${sample_id}.clean.R1_bismark_bt2_pe.deduplicated.CX_report.txt \
  -o ${sample_id}.clean.R1_bismark_bt2_pe.deduplicated.CX_report.CGmap.gz

  $software/cgmaptools mstat \
  -i ${sample_id}.clean.R1_bismark_bt2_pe.deduplicated.CX_report.CGmap.gz \
  -c ${ncpus} -f pdf > ${sample_id}.clean.R1_bismark_bt2_pe.deduplicated.CX_report.CGmap_mstat.data

  $software/cgmaptools mtr \
  -i ${sample_id}.clean.R1_bismark_bt2_pe.deduplicated.CX_report.CGmap.gz \
  -r $REFERENCE_DIR/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.111.sort.bed \
  -o ${sample_id}.clean.R1_bismark_bt2_pe.deduplicated.CX_report.CGmap_mtr.gz

  gzip CHG_context_${sample_id}.clean.R1_bismark_bt2_pe.deduplicated.txt
  gzip CHH_context_${sample_id}.clean.R1_bismark_bt2_pe.deduplicated.txt
  gzip CpG_context_${sample_id}.clean.R1_bismark_bt2_pe.deduplicated.txt
  gzip ${sample_id}.clean.R1_bismark_bt2_pe.deduplicated.CX_report.txt

elif [ "$paired_end" -eq 1 ]; then
  bam_file="${out_dir}/${sample_id}/bismark/${sample_id}.clean_bismark_bt2.deduplicated.bam"
  CX_file="${sample_id}.clean_bismark_bt2.deduplicated.CX_report.txt"
  echo "Input bam files: $bam_file"
  echo "Input CX files: $CX_file"

  bam2nuc \
  --genome_folder $ref_dir/bismark.index \
  --dir ${out_dir}/${sample_id}/bismark \
  --samtools_path /p300s/zhaowm_group/tangbx/software/miniconda3/envs/RNAseq_E4/bin/ \
  ${out_dir}/${sample_id}/bismark/${sample_id}.clean_bismark_bt2.deduplicated.bam

  $software/cgmaptools convert bismark2cgmap \
  -i ${sample_id}.clean_bismark_bt2.deduplicated.CX_report.txt \
  -o ${sample_id}.clean_bismark_bt2.deduplicated.CX_report.CGmap.gz

  $software/cgmaptools mstat \
  -i ${sample_id}.clean_bismark_bt2.deduplicated.CX_report.CGmap.gz \
  -c ${ncpus} -f pdf > ${sample_id}.clean_bismark_bt2.deduplicated.CX_report.CGmap_mstat.data

  $software/cgmaptools mtr \
  -i ${sample_id}.clean_bismark_bt2.deduplicated.CX_report.CGmap.gz \
  -r ${ref_dir}/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.111.sort.bed \
  -o ${sample_id}.clean_bismark_bt2.deduplicated.CX_report.CGmap_mtr.gz

  gzip CHG_context_${sample_id}.clean_bismark_bt2.deduplicated.txt
  gzip CHH_context_${sample_id}.clean_bismark_bt2.deduplicated.txt
  gzip CpG_context_${sample_id}.clean_bismark_bt2.deduplicated.txt
  gzip ${sample_id}.clean_bismark_bt2.deduplicated.CX_report.txt

else
    echo "ERROR: $bam_file or $CX_file does not exist!"
    exit 1
fi

echo "All programs completed!"
