#!/bin/bash -e

# @File       :bulkRNA-seq_E4.sh
# @Time       :2024/5/6 16:12
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :bulkRNA_pipeline
# @Version    :1.0.0
# @Description:measure RNA abundance
# @Usage      :bash bulkRNA-seq_E4.sh <fastq_dir> <out_dir> <sample_id> <ncpus> <ramGB>

###software vision###
#fastp:0.23.2
#kallisto:0.46.0
#STAR:2.73a
#rsem:1.3.3

###parameter###
ref_dir=/path/reference/UU_Cfam_GSD_1.0
fastq_dir=$1
out_dir=$2
sample_id=$3
ncpus=$4
ramGB=$5
fwd_prob=0.5
#############

# Check for the existence of files ending with '_1.fastq.gz' or '_2.fastq.gz'
PE_raw=$(find "${fastq_dir}/${sample_id}/reads/" -type f \( -name "${sample_id}_1.fastq.gz" -o -name "${sample_id}_2.fastq.gz" \))
# Check for the existence of files ending with '.fastq.gz'
SE_raw=$(find "${fastq_dir}/${sample_id}/reads/" -type f -name "${sample_id}.fastq.gz")

###run fastp###
mkdir -p ${out_dir}/${sample_id}/fastp
###Determine sequencing strategy and execute code
if [[ -n "$PE_raw" && -z "$SE_raw" ]]; then
    # clean data for PE reads
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

elif [[ -z "$PE_raw" && -n "$SE_raw" ]]; then
    # clean data for SE reads
    echo "...run fastp SE pipeline..."
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


# Check for the existence of files ending with 'clean.R1.fastq.gz' or 'clean.R2.fastq.gz'
PE_clean=$(find "${out_dir}/${sample_id}/fastp/" -type f \( -name "${sample_id}.clean.R1.fastq.gz" -o -name "${sample_id}.clean.R2.fastq.gz" \))
# Check for the existence of files ending with 'clean.fastq.gz'
SE_clean=$(find "${out_dir}/${sample_id}/fastp/" -type f -name "${sample_id}.clean.fastq.gz")

### run kallisto###
mkdir -p ${out_dir}/${sample_id}/kallisto_gene
mkdir -p ${out_dir}/${sample_id}/kallisto_transcript
mkdir -p ${out_dir}/${sample_id}/star

if [[ -n "$PE_clean" && -z "$SE_clean" && "$paired_end" -eq 0 ]]; then
    echo "Input files: $PE_clean"
    kallisto quant \
    -t ${ncpus} --fusion --plaintext \
    -i ${ref_dir}/kallisto.index/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.111.gene.fa.idx \
    -o ${out_dir}/${sample_id}/kallisto_gene \
    ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R1.fastq.gz \
    ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R2.fastq.gz

    kallisto quant \
    -t ${ncpus} --fusion --plaintext \
    -i ${ref_dir}/kallisto.index/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.111.transcript.fa.idx \
    -o ${out_dir}/${sample_id}/kallisto_transcript \
    ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R1.fastq.gz \
    ${out_dir}/${sample_id}/fastp/${sample_id}.clean.R2.fastq.gz

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

elif [[ -z "$PE_clean" && -n "$SE_clean" && "$paired_end" -eq 1 ]]; then
    echo "Input files: $SE_clean"
    kallisto quant \
    -t ${ncpus} --fusion --plaintext --single -l 200 -s 20 \
    -i ${ref_dir}/kallisto.index/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.111.gene.fa.idx \
    -o ${out_dir}/${sample_id}/kallisto_gene \
    ${out_dir}/${sample_id}/fastp/${sample_id}.clean.fastq.gz

    kallisto quant \
    -t ${ncpus} --fusion --plaintext --single -l 200 -s 20 \
    -i ${ref_dir}/kallisto.index/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.111.transcript.fa.idx \
    -o ${out_dir}/${sample_id}/kallisto_transcript \
    ${out_dir}/${sample_id}/fastp/${sample_id}.clean.fastq.gz

    STAR \
    --genomeDir ${ref_dir}/star.index \
    --readFilesIn ${out_dir}/${sample_id}/fastp/${sample_id}.clean.fastq.gz \
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
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate  \
    --quantMode TranscriptomeSAM \
    --sjdbScore 1 \
    --limitBAMsortRAM ${ramGB}000000000

else
    echo "ERROR: Clean sequencing files does not exist or has a wrong name!"
    exit 1
fi

###run rsem###
mkdir -p ${out_dir}/${sample_id}/rsem
bam_file="${out_dir}/${sample_id}/star/${sample_id}.gsd.Aligned.sortedByCoord.out.bam"

# Checking for the BAM files
if [ -f "$bam_file" ]; then
    STAR \
    --runMode inputAlignmentsFromBAM \
    --inputBAMfile ${out_dir}/${sample_id}/star/${sample_id}.gsd.Aligned.sortedByCoord.out.bam \
    --outWigType bedGraph \
    --outWigStrand Stranded \
    --outWigReferencesPrefix chr \
    --outFileNamePrefix ${out_dir}/${sample_id}/star/${sample_id}.gsd.

    rsem-calculate-expression --version

    if [ "$paired_end" -eq 0 ]; then
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

    elif [ "$paired_end" -eq 1 ]; then
        rsem-calculate-expression \
        --bam \
        --estimate-rspd \
        --seed 12345 \
        -p ${ncpus} \
        --no-bam-output \
        --forward-prob ${fwd_prob} \
        ${out_dir}/${sample_id}/star/${sample_id}.gsd.Aligned.toTranscriptome.out.bam \
        $ref_dir/rsem.index/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa \
        ${out_dir}/${sample_id}/rsem/${sample_id}_rsem

    else
        echo "ERROR：Please input PE or SE sequencing files."
        exit 1
    fi
else
    echo "ERROR: $bam_file does not exist!"
    exit 1
fi

echo "All programs completed!"

