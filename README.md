# Current Version: L1ME-AID (v1.4.0-beta)
***If you are looking for an earlier version (e.g., v1.0.0-beta) check the previousVersions folder.*** </br>
-Note: Starting at Version 1.3.0-beta multiple column names have been adjusted. This is to better reflect that some sequences are not TEs. </br>
-Note: Beginning at Version 1.1.0-beta L1ME-AID will check for TSD sequence if you provide the reference genome as -g and name your sequences as 'chromosome-position-anythingElse' (e.g., chr1-1002321-whatever). TSD check will only happen if you give -g a path to a reference file otherwise this functionality is skipped. PLEASE be aware this is only tested for TSD finding for INSERTIONS only. It will run on all rows as L1ME-AID expects only insertions but if you are just using L1ME-AID to check DEL calls this can be deceiving if you try and run TSDs as well. <br>

***How can I cite?***</br>
We are currently working on polishing L1ME-AID and on a manuscript. Currently, if you are using L1ME-AID and wanting to cite it please use Logsdon, G.A., Ebert, P., Audano, P.A. et al. Complex genetic variation in nearly complete human genomes. Nature (2025). https://doi.org/10.1038/s41586-025-09140-6<br>

## L1ME-AID (L1 Mediated Annotation and Insertion Detector)
<p align="center">
<img src="limeaid.jpeg?raw=true" width="350" height="350">
</p>
L1ME-AID is a simple python script that takes two files as input:<br>
  1) Fasta file of structural variation sequences you want to check for L1 mediated transposition events (Mobile element insertions <i>Alu</i> elements, L1, SVA, etc.)<br>
  2) RepeatMasker .out file that was produced from running RepeatMasker on the Fasta file.<br>

# Installation
Download the code and run the python file (see HOW TO RUN below).

# Requirements:
This program was tested on Linux (Mint 21.2). It has not been tested on Mac/PC.

## Hardware:
A standard computer with enough RAM to support the in-memory operations/dataset you provide.

## Python Dependencies 
(versions are listed below, but the functions utilized are very basic it will most likely work with most versions. Additionally, many of these are built-in to python (e.g., argparse):<br>
  1) pandas - v2.2.1<br>
  2) pysam -  v0.22.0<br>
  3) os - Built-in (os)<br>
  4) ast - Built-in (ast)<br>
  5) numpy - v1.26.4<br>
  6) json - v2.0.9<br>
  7) collections - Built-in (collections)<br>
  8) from Bio.Seq import Seq - v1.83<br>
  9) from Bio import SeqIO -v1.83<br>
  10) more_itertools - v10.2.0 <- you can get more_itertools through pip if Anaconda is no longer available to you<br>
  11) argparse - v1.1<br>
  12) from functools import reduce<br>

# HOW TO RUN:

## Minimum Code Necessary:
`python limeaid.py -i /your/path/to/fastafile/demo.fasta -r /your/path/to/repeatmaskerfile/demo.fasta.out -o /your/path/to/outfile/demoOutFile.tsv`
- Provide an input Fasta file (-i), the Repeatmasker .out file (-r), and the path to the tsv that will be exported (-o).

## Demo:
Check out the demo files for a quick test of L1ME-AID:</br>

You can see if it works for you by running this code: `python limeaid.py -i demo.fasta -r demo.fasta.out -o ./demo.outFile`</br>
The expected run time is less than 1 second for this demo. You should expect to get the same file as the demoOutFile.tsv provided. 

## Other options

  a) -rd, Repeat Divergence, default=20% divergence; 20% maximum divergence for reading through all the repeatmasker annotations (I would leave this but you can change it by providing an integer (-rd 10 - for 10% divergence maximum)). <br>
  b) -d, divergences, default=6,15,15 ; 6% for <i>Alu</i> elements, 15% for L1s and 15% for SVAs (provide a each number in this list format x,y,z, EXAMPLE: -d 5,10,15 must have three numbers and each number corresponds to Alu,L1,SVA)<br>
  c) -l, length cutoffs, default=500,10000,10000; 500bp for <i>Alu</i> elements, 10kbp for L1, 10kbp for SVA (provide a each length in this list format x,y,z, EXAMPLE: -l 300,8000,3000) <br>
  d) -u, upper length cutoff, default=50000, if a sequence is longer than this don't even bother<br>
  e) -t, maximum tail start position, default=50; this will allow the start of an elements tail to begin at position 50 in the sequence (feel free to change as desired)<br>
  f) -yl, allowable LINE subfamilies, default='L1HS,L1P1,L1PA1,L1PA2'; These are the L1 subfamilies that L1ME-AID will consider to be active. If for example an L1PA3 is annotated this will be noted as OLDER LINE SUBFAMILY<br>
  g) -g, path to the reference genome utilized during SV calling (e.g., '/home/yourname/myfolder/hg38.fa'). *Note: Only provide a reference genome here if you want to have L1ME-AID try to call TSDs. If you do this you need to name the sequences in your fasta in a specific format (chromosome-position-anything else: chr1-101101-mysequence). L1ME-AID uses this chromosome and position information to look for the coordinate sequence that matches in your SV sequence. 

## Pay Attention
Please pay attention to the following columns:<br>
  a) ID = the IDs need to be the same in the fasta file as the repeatmasker out file<br>
  b) Element_Percentage: What proportion of the sequence was annotated as specific elements<br>
  c) Element_Designation: The type of element that makes up the largest proportion of the sequence<br>
  d) FILTER_RESULTS: This will provide you the results from the filters (divergences, lengths, etc.). It will say 'Good_Row' if the element passes the checks. HOWEVER, it will always say 'Good_Row' for sequences that are not <i>Alu</i> elements, L1s, and SVAs. Additionally, this DOES NOT take into account other columns (Tail_Type, Tail_Length, Unique_Element_Count). If you are looking for YOUNG mobile element insertions then this column is important. If you are not or are more flexible on what you want to call feel free to just IGNORE THIS COLUMN. <br>
  e) Tail_Begins: Where in the sequence the tail begins (this number corresponds to a T tail so we are flipping elements and numbering from the 5' side - just in case you were wondering why your A tail is at position 1).<br>
  f) Tail_Type: This will try and tell you if there is a tail possibly what type is present. If you are looking for active mobile elements there should be a tail. FILTER_RESULTS will say 'Good_Row' even if there is NO TAIL. Make sure to filter your data by both columns.<br>
  g) Tail_Length: Tries to estimate the length of the tail. Finding tails uses a 5bp seed so the minimum should be 5bp.<br>
  h) Tail_Seed_Hits: How many unique 5bp seeds were found (not too important for the user).<br>
  i) Unique_Element_Count: This will tell you if multiple subfamilies were identified in a sequence (Ex: <i>AluS</i>, <i>AluY</i>). This is more important for <i>Alu</i> elements. The code doesn't care about this column but it's important for the user to decide. Usually it's not a big deal if L1s or SVAs have multiple subfamilies as RepeatMasker will sometimes annotate pieces of the element separately. <br>

# Please Remember:
This program will provide you all of the output no matter the results. You can always choose to ignore the results if you think something is wrongly annotated. It is only here to help you in finding active elements mobilized by L1 machinery. This is still in beta so it might not be perfect. I usually perform manual curation of calls after. 

# License
This project is covered under the GNU Lesser General Public License, version 3.0 (LGPL-3.0)

# Log
September 12, 2025: After accidentally setting L1ME-AID up on a compute cluster that was pointing to a python2 installation I found some bugs in element proportion calculations. This shouldnt break if installed with python3 but made some quick changes to help ensure if others accidentally use python2 it shouldnt break. 
September 12, 2025: reuploaded 1.3.4. On a compute cluster a comment was acting as a non-ASCII character so I removed that comment. All functions stayed the same. I also removed the need for tqdm. Trying to start limiting external packages required to make things simpler. 
August 19, 2025: v1.3.4-beta uploaded with some small changes. The Filter column was renamed Flags, and 'Good_Row' value within that column is now 'No_Flags'. This is to better reinforce these are just suggestions to check not 100% rely on. <br> 
August 13, 2025: Complete rework of the element divergence estimates. <br>
July 1, 2025: Sometimes repeatmasker does not include an ID column (#14). Added a fix to work around this. <br>
June 20, 2025: Fixed minor bug where TWIN PRIMING elements orientation was not caught with Tail Checks. <br>
June 15, 2025: Major updates to multiple functions. For example, Divergences now are length-adjusted where two+ of the same element annotations are noted in a sequence, orientations should be more accurate.<br>
May 29, 2025: Fixed minor bug. <br>
March 4, 2025: L1ME-AID v1.2.0-beta was uploaded. Updated RepeatMaskerPatternFinder function as well as the iterative function application at the end. <br>
January 17, 2025: L1ME-AID exports the file as a TSV instead of CSV now. <br>
December 09, 2024: Fixed a bug where non-MEIs without TSDs were breaking the code. Version 1.1.1 was uploaded. Seems to work. <br>
November 19, 2024: Verson 1.1.0-beta uploaded, added the -g functionality to call TSDs that are present within the SV sequence and the reference genome. 
