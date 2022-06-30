#!/bin/bash

#SBATCH -n 10
#SBATCH --mem 10G

project_path=$1
$basespace_path=$2
run=$3
sample=$4

echo Started merging fastq files for ${sample}.

cat $basespace_path/Projects/${run}/Samples/${sample}/Files/*_R1_001.fastq.gz > $project_path/fastq/${sample}_R1.fastq.gz
cat $basespace_path/Projects/${run}/Samples/${sample}/Files/*_R2_001.fastq.gz > $project_path/fastq/${sample}_R2.fastq.gz


echo Finished merging fastq files for ${sample}.

echo Started STAR alignment for ${sample}.

STAR --runThreadN 8 --genomeDir /storage/westbrook/genomes/hg38/sequence/star_index_gencode_v36 \
--readFilesIn ''${project_path}'/fastq/'${sample}'_R1_merged.fastq.gz' ''${project_path}'/fastq/'${sample}'_R2_merged.fastq.gz' \
--outFileNamePrefix ''${project_path}'/bam/'${sample}'_' --readFilesCommand zcat \
--sjdbGTFfile /storage/westbrook/genomes/hg38/annotation/gencode_v36_annotation.gtf  \
--outReadsUnmapped Fastx --alignSJDBoverhangMin 1 --alignSJoverhangMin 8 \
--outSAMtype BAM SortedByCoordinate --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
--outFilterMultimapNmax 100 --winAnchorMultimapNmax 100

echo Finished STAR alignment for ${sample}.

echo Started TEcount quantification for ${sample}.

TEcount -b ${project_path}/bam/${sample}_Aligned.sortedByCoord.out.bam --sortByPos \
--stranded reverse --GTF /storage/westbrook/genomes/hg38/annotation/hg38_GENCODE_v36.gtf \
--TE /storage/westbrook/genomes/hg38/annotation/GRCh38_GENCODE_rmsk_TE.gtf \
--project ${project_path}/TEcount/${sample}

echo Finished TEcount quantification for ${sample}.
