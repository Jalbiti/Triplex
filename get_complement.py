#!/usr/bin/env python

import argparse
import sys
import pandas as pd
import numpy as np

inFile = sys.argv[1]
outFile = sys.argv[2]

# Strictly Parallel Triplex Triads
def compl_seq(seq):
        complementary = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
        return ''.join([complementary[x] for x in seq])

triplex = pd.read_csv(inFile, sep ="\t", header=None)
seqs = triplex[0]

compl = [compl_seq(seq) for seq in seqs]

fout = open(outFile, "w")
c = 0
for seq in compl:
  fout.write(seq + "\n")
  c += 1

fout.close()


