#!/usr/bin/env python3
#  GENERAL USAGE: python ./scripts/Divide_Reference.fa_by_Contig.py ./data/reference/reference.fa

import os
from os import path
import sys
import argparse

parser = argparse.ArgumentParser("Script split a combined fasta file with multiple sequences into one fasta file per sequence.  Output files are named after the sequence.")
parser.add_argument("input", help="The multi-sequence fasta file to split")
parser.add_argument("-o", "--outdir", required=False, default="split_out", help="Output directory for split fasta files")
args = parser.parse_args()

infile=open(args.input)
outdir=args.outdir

if not os.path.exists(outdir):
    os.makedirs(outdir)

num_contigs = 'grep "^>" %s  | wc -l' % (args.input)
print("The number of contigs in reference fasta:")
os.system(num_contigs)


opened = False # Assume outfile is not open

for line_ref in infile:
    if line_ref[0] == ">": # If line begins with ">"
        if(opened): 
            outfile.close() # Will close the outfile if it is open (see below and follow loop)
        opened = True # Set opened to True to represent an opened outfile
        contig_name = line_ref[1:].rstrip() #Extract contig name: remove ">", extract contig string, remove any spaces or new lins following file
        print("contig: " + contig_name)
        outfilepath=os.path.join(outdir, str(contig_name) + ".fa")
        print("Name of outfile: " +  outfilepath)
        outfile=open(outfilepath, 'w')
    outfile.write(line_ref) # Write the line to the file. If line started with ">" a new output file was created and ready to be written to.
outfile.close()

