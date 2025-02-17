#!/bin/bash
#SBATCH -p AMDMEM
#SBATCH -o feats.o
#SBATCH -e feats.e

module load ANACONDA2
conda activate triplex

inFile=$1
#k=$2
#chr=$3
chr=$2

#Rscript features.R $inFile $k $chr 
Rscript features.R $inFile $chr

