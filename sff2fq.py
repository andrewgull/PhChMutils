#!/usr/bin/env python3

from Bio import SeqIO
import sys

sff_file = sys.argv[1]
fq_file = sys.argv[2]

sff = [r for r in SeqIO.parse(sff_file, 'sff')]

SeqIO.write(sff, fq_file, 'fastq')