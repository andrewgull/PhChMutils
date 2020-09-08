#!/home/andrei/miniconda3/bin/python3

from Bio import AlignIO
import sys

infile = sys.argv[1]

def get_extension_length(word):
    if word[-5:] == 'fasta':
        ext = 5
    elif word[-2:] == 'fa':
        ext = 2
    elif word[-3:] == 'fas':
        ext = 3
    else:
        ext = 0
    return ext

def write_aln(filename, extlen):
    if extlen > 0:
        outname = infile[:-n]+'phy'
        count = AlignIO.convert(infile, 'fasta', outname, 'phylip-relaxed')
    else:
        outname = infile+'phy'
        count = AlignIO.convert(infile, 'fasta', outname, 'phylip-relaxed')
    return (count, outname)

n = get_extension_length(infile)
outfile = write_aln(infile, n)[1]

print('Done ----> %s' %outfile)

