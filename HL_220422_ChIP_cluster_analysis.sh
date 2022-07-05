#!/usr/bin/bash

path=$1
sample=$2

## get gene list from table
#awk -F '\t' '{if($2=="${sample}") {print}}' $path/Timex_foldchange_clusters.txt | cut -f1 > $path/cluster_${sample}.gene.txt

grep -f $path/cluster_${sample}.gene.txt /storage/westbrook/genomes/hg38/annotation/gencode_v36_annotation.gtf | awk '!a[$1 $4]++' > $path/cluster_${sample}.gtf

computeMatrix reference-point -S \
  ~/side_analysis/timex_chip/bigwig/MYC_4H_DMSO_1.bw \
  ~/side_analysis/timex_chip/bigwig/MYC_4H_DMSO_2.bw \
  ~/side_analysis/timex_chip/bigwig/MYC_4H_TAM_1.bw \
  ~/side_analysis/timex_chip/bigwig/MYC_4H_TAM_2.bw \
  ~/side_analysis/timex_chip/bigwig/MYC_12H_DMSO_1.bw \
  ~/side_analysis/timex_chip/bigwig/MYC_12H_DMSO_2.bw \
  ~/side_analysis/timex_chip/bigwig/MYC_12H_TAM_1.bw \
  ~/side_analysis/timex_chip/bigwig/MYC_12H_TAM_2.bw \
  -R $path/cluster_${sample}.gtf \
  --referencePoint TSS -a 3000 -b 3000 \
  -o $path/cluster_${sample}

plotHeatmap -m $path/cluster_${sample} \
  -o $path/cluster_${sample}.pdf \
  --sortUsingSamples 1 --heatmapHeight 15 --colorMap hot \
  --whatToShow 'plot, heatmap and colorbar' \
  --samplesLabel 4H_DMSO_1 4H_DMSO_2 4H_TAM_1 4H_TAM_2 12H_DMSO_1 12H_DMSO_2 12H_TAM_1 MYC_12H_TAM_2
