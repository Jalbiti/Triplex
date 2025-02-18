# Triplex


This repository contains the code used in the article "Systematic study of hybrid triplex topology and stability suggests a general triplex-mediated regulatory mechanism" by Genna et al. published in NAR, to study potential TFO in various datasets.

### Bioinformatics pipeline
**The pipeline** consists on 4 steps:
-  Finding candidate sequences (e.g. `find_seqs.py`)
-  Extracting the candidate target sites (e.g. `get_complement.py`)
-  Aligning the triplex target sites (TTS) to a reference genome (e.g. `align.sh`)
-  Comparing our mapped sites to other elements, such as nucleosome positioning (e.g. `features.sh` and `nuc_plots.R`)

An exemplary pipeline can be found in `pipeline.sh` to find candidate triplexes from annotated miRNAs and lncRNAs.
