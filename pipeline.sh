#!/bin/bash

python find_seqs.py

python get_complement.py miRNA.out miRNA_c.out
python get_complement.py lncRNA.out lncRNA_c.out

# Now align complements

#As string search takes too long we will just align with STAR
python seqsToFasta.py -i miRNA_c.out -o miRNA_TTS.fa
python seqsToFasta.py -i lncRNA_c.out -o lncRNA_TTS.fa

sbatch align.sh miRNA_TTS.fa miRNA_TTS
sbatch align.sh lncRNA_TTS.fa lncRNA_TTS

# Get features with features.sh and nucleosome plots
#sbatch features.sh
#Rscript nuc_plots.R $1
