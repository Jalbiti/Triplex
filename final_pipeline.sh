#!/bin/bash

#SBATCH -J pipe
#SBATCH -e pipe.e
#SBATCH -o pipe.o
#SBATCH --partition=XESH
#SBATCH --mem=75GB 
#SBATCH --time=03:00:00
#SBATCH --ntasks=10
##SBATCH --exclusive
##SBATCH -p MPI
##SBATCH --cpus-per-task=8


module load STAR/2.7.5a R samtools GCC

fasta=$1

#########STAR Alignment to human genome
STAR --genomeDir /orozco/projects/NucleosomeDynamics/Alba/NCBI/human/STAR/ \
 --runThreadN 10 \
 --readFilesIn $fasta \
 --outFileNamePrefix outSTAR_ \
 --outFilterMultimapNmax 20 \
 --quantMode TranscriptomeSAM GeneCounts \
 --outSAMtype BAM SortedByCoordinate \
 --outSAMattributes All

# Create BAM index
samtools index outSTAR_Aligned.sortedByCoord.out.bam
samtools fasta outSTAR_Aligned.sortedByCoord.out.bam > outSTAR_.fa


#########Triplex Stability Predictor
module load ANACONDA2
conda activate triplex

python FindTFO_Test1.py  -i outSTAR_.fa  -t DNA -o outputSTAR.fa | grep OK > results.out

#########BWA Alignment of Shorter Fragments

# Align Fasta (after applying predictor)
/orozco/projects/NucleosomeDynamics/Alba/NCBI/BWA/bwa-0.7.17/bwa aln /orozco/projects/NucleosomeDynamics/Alba/NCBI/human/BWA/bwa_db outputSTAR.fa > first.bwa
echo "Done with BWA Part1"

# Align complemet seqs of fragments to analyze mapping to other regions
/orozco/projects/NucleosomeDynamics/Alba/NCBI/BWA/bwa-0.7.17/bwa samse /orozco/projects/NucleosomeDynamics/Alba/NCBI/human/BWA/bwa_db first.bwa outputSTAR.fa > first.se.sam
echo "Done with BWA Part2"

# Generate BAM files
samtools view -S -b first.se.sam > first.se.bam
echo "Done with Part3"

# Sort BAM File
samtools sort first.se.bam -o first.se.sorted.bam
echo "Done with Part4"

# Index BAM file
samtools index first.se.sorted.bam
echo "Done with Part5"

# Count mapped reads (using Rscript)
Rscript fCounts.R first.se.sorted.bam > Multi_filter
Rscript fCounts2.R first.se.sorted.bam > noMulti_filter

############## PRE FILTER
#Rscript fCounts.R /orozco/projects/NucleosomeDynamics/Alba/NCBI/STAR/results/org_align/EMTAB8300Aligned.sortedByCoord.out.bam > Multi_org
#Rscript fCounts2.R /orozco/projects/NucleosomeDynamics/Alba/NCBI/STAR/results/org_align/EMTAB8300Aligned.sortedByCoord.out.bam > noMulti_org

############## POST FILTER
#Rscript fCounts.R ../../human/BWA/EMTAB8300_C_multi.bam > Multi_filt
#Rscript fCounts2.R ../../human/BWA/EMTAB8300_C_multi.bam > noMulti_filt



