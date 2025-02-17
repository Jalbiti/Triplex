# -*- coding: utf-8 -*-

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="\tSeqs to Fasta...\n",\
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', required=True, metavar='input_file',type=str, help="Input sequences")
parser.add_argument('-o', required=True, metavar='output_file',type=str, help="Output sequences in FASTA format")

args = parser.parse_args()

inFile = args.i
outFile = args.o
records = pd.read_csv(inFile, header = None)
sequences = []
for rec in range(len(records)):
    sequences.append(records.iloc[rec][0])

fout = open(outFile, 'w')
with open(inFile) as fp:
    cnt = 1
    for seq in sequences:
        fout.write(">Seq {}\n".format(cnt))
        fout.write(str(seq) + "\n")
        #fout.write("\n")
        cnt += 1
   
fout.close()

