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

