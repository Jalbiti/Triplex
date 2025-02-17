#!/bin/bash

#########Triplex Stability Predictor

fasta=$1
out_name=$2

#########BWA Alignment of Shorter Fragments

STAR --genomeDir STAR/hg38_new \
 --runThreadN 10 \
 --readFilesIn $fasta \
 --outFileNamePrefix $out_name. \
 --outFilterMultimapNmax 20000000000 \
 --winAnchorMultimapNmax 20000000000 \
 --alignTranscriptsPerReadNmax 100000 \
 --limitBAMsortRAM 72753112415 \
 --quantMode TranscriptomeSAM \
 --outSAMtype BAM SortedByCoordinate \
 --outSAMattributes All

# --limitOutSJcollapsed 100000000 \

# Create BAM index and extract mapping per chr
mv $out_name.Aligned.sortedByCoord.out.bam temp_bam.bam
samtools index temp_bam.bam
samtools view -b temp_bam.bam chr1 > $out_name.chr1.bam
samtools view -b temp_bam.bam chr2 > $out_name.chr2.bam
samtools view -b temp_bam.bam chr3 > $out_name.chr3.bam
samtools view -b temp_bam.bam chr4 > $out_name.chr4.bam
samtools view -b temp_bam.bam chr5 > $out_name.chr5.bam
samtools view -b temp_bam.bam chr6 > $out_name.chr6.bam
samtools view -b temp_bam.bam chr7 > $out_name.chr7.bam
samtools view -b temp_bam.bam chr8 > $out_name.chr8.bam
samtools view -b temp_bam.bam chr9 > $out_name.chr9.bam
samtools view -b temp_bam.bam chr10 > $out_name.chr10.bam
samtools view -b temp_bam.bam chr11 > $out_name.chr11.bam
samtools view -b temp_bam.bam chr12 > $out_name.chr12.bam
samtools view -b temp_bam.bam chr13 > $out_name.chr13.bam
samtools view -b temp_bam.bam chr14 > $out_name.chr14.bam
samtools view -b temp_bam.bam chr15 > $out_name.chr15.bam
samtools view -b temp_bam.bam chr16 > $out_name.chr16.bam
samtools view -b temp_bam.bam chr17 > $out_name.chr17.bam
samtools view -b temp_bam.bam chr18 > $out_name.chr18.bam
samtools view -b temp_bam.bam chr19 > $out_name.chr19.bam
samtools view -b temp_bam.bam chr20 > $out_name.chr20.bam
samtools view -b temp_bam.bam chr21 > $out_name.chr21.bam
samtools view -b temp_bam.bam chr22 > $out_name.chr22.bam

echo "DONE"

