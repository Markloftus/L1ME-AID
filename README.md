# L1ME-AID (beta version)
## L1ME-AID (L1 Mediated Annotation and Insertion Detection)

L1ME-AID is a simple python script that takes two files as input:
  1) Fasta file of structural variation sequences you want to check for L1 mediated transposition events (Mobile element insertions Alus, L1, SVA, etc.)
  2) RepeatMasker .out file that was produced from running RepeatMasker on the Fasta file

## You will need to be able to import these libraries:
  1) import matplotlib.pyplot as plt
  2) import pandas as pd
  3) import pysam
  4) import os
  5) import ast
  6) import numpy as np
  7) import json
  8) import collections
  9) from tqdm import tqdm
  10) from Bio.Seq import Seq
  11) from Bio import SeqIO
  12) import more_itertools as mit
  13) import argparse as ap

# HOW TO RUN:

## Minimum Code Necessary:
python limeaid.py -i /home/mark/Desktop/MEI_Group/HGSVC2/insertions/hgsvc2_INS.fasta -r /home/mark/Desktop/MEI_Group/HGSVC2/insertions/repeatmasker/hgsvc2_INS.fasta.out -o /home/mark/Desktop/test.csv
- Provide an input Fasta file (-i), the Repeatmasker out file (-r), and the path to the csv that will be exported (-o).

## Other options

  a) -rd, Repeat Divergence, default=20% divergence; 20% maximum divergence for reading through all the repeatmasker annotations (I would leave this but you can change it by providing an integer.
  b) -d, divergences, default=6,15,15 ; 6% for Alus, 15% for L1s and 15% for SVAs (provide a each number in this list format x,y,z)
  c) -l, length cutoffs, default=500,10000,10000; 500bp for Alus, 10kbp for L1, 10kbp for SVA (provide a each length in this list format x,y,z)
  d) -u, upper length cutoff, default=50000, if a sequence is longer than this don't even bother
  e) -t, maximum tail start position, default=50; this will allow the start of an elements tail to begin at position 50 in the sequence (feel free to change as desired)
  f) -yl, allowable LINE subfamilies, default=L1HS,L1PA2; These are the L1 subfamilies that L1ME-AID will consider to be active. If for example an L1PA3 is annotated this will be noted as OLDER LINE SUBFAMILY

# Remember:
This program will provide you all of the output no matter what the results are. You can always choose ignore the results if you think something is wrongly annotated. It is only here to help you in finding active elements mobilized by L1 machinery.

## Pay Attention
Please pay attention to the following columns:
  a) 
