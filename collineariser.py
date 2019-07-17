#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
collinearise aligment made by MAUVE
MAUVE makes merged alignment file containing each LCB for each genome 
that actually contains it (e.g. length of LCB blocks can be smaller than number of genomes)
The idea is to concatenate LCB blocks in a horizontal way

Created on Mon Jun 10 18:32:45 2019

@author: andrei guliaev

"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import sys
import os

help_message = "Collinearise aligment produced by MAUVE\n" \
               "MAUVE makes merged alignment file containing each LCB for each genome\n" \
               "that actually contains it (e.g. length of LCB blocks can be smaller than number of genomes)\n" \
               "The idea is to concatenate LCB blocks in a horizontal way.\n" \
               "NB: sequence ID must not contain 'fasta' word (except if it is extension)" \
               "Arguments:\n" \
               "1. - alignment file from MAUVE (usually has no extension)\n" \
               "2. - number of genomes\n" \
               "3. - output file name."

if len(sys.argv) < 2:
    print(help_message)
    sys.exit()

# filename = 'out_clean.fasta' <- alignment file which is produced by MAUVE 
# has no extension
# by default MAUVE out file contains commented lines
# which must be removed in advance

filename = sys.argv[1]
# ngenomes = 41
ngenomes = int(sys.argv[2])
outname = sys.argv[3]

# remove commented lines
os.system("sed '/#/d' %s > %s_aln.fasta" % (filename, filename))
filename += '_aln.fasta'
os.system("sed -i 's/=//g' %s" % filename)

# mauve = [rec for rec in SeqIO.to_dict(SeqIO.parse(filename,'fasta'))]
# read fasta as dict
mauve = SeqIO.to_dict(SeqIO.parse(filename,'fasta'))
# each record in this quasi-alignment has its ID/order number in the beggining,
# like '1:116500-116934'

x = len(mauve[list(mauve.keys())[0]]) # alignment length

def get_first_number(string):
    number = int(re.findall(pattern="[0-9]*[:,_]", string=string)[0][:-1])
    return number

# GET BORDER INDEXES
print('1. GET BORDER INDEXES')
# v - list of the dict keys
#v = ['1:_','2:_','3:_','4:_','1:-','2:-','3:-','2:+','4:+','2:*','3:*','4:*']
#s = ['att', 'att', 'attt', 'attt', 'acc', 'acc', 'acc', 'ggg', 'ggg', 'gga', 'gga', 'gga']

#vseq = [SeqRecord(Seq(item[1]), id=item[0]) for item in zip(v,s)]
#vd = SeqIO.to_dict(vseq)
vd = mauve
v = list(vd.keys())

border_idx = list()
n = 1
#print("LCB %i - %i" %(n, get_first_number(v[0])))
sub = [v[0]]

for i in range(1, len(v)):
    if get_first_number(v[i]) < get_first_number(v[i-1]) or len(vd[v[i]]) != len(vd[v[i-1]]):
        border_idx.append(i)
        n += 1

# =============================================================================
# for i in range(1, len(v)):
#     if get_first_number(v[i]) > get_first_number(v[i-1]):
#         print("LCB %i - %i" %(n, get_first_number(v[i])))
#     else:
#         border_idx.append(i)
#         print("index in v - %i" %(i))
#         n += 1
#         print("LCB %i - %i" %(n, get_first_number(v[i])))
# print("index in v - %i" %len(v))
# =============================================================================

border_idx.append(len(v))
        
# SLICE LIST BY INDEXES
print('2. SLICE LIST BY INDEXES')
border_idx = [0] + border_idx

# GET INDIVIDUAL LCBs
print('3. GET INDIVIDUAL LCBs')
sliced_lst = [v[slice(border_idx[i], border_idx[i+1])] for i in range(len(border_idx)-1)]

# INSERT EMPTY RECORDS
print('4. INSERT EMPTY RECORDS')
lacking_rec_all = list()
for item in sliced_lst:
    existing = [get_first_number(subitem) for subitem in item]
    lacking = set(range(1,ngenomes+1)) - set(existing)
    #print(lacking)
    lacking_rec = [str(lack) + '_empty' for lack in lacking]
    item += lacking_rec
    lacking_rec_all += lacking_rec
    item.sort()

# all IDs of lacking recs are stored here
lacking_rec_all = set(lacking_rec_all)


# COLLECT SEQRECORDS according to sliced list of IDs
print('5. COLLECT SEQRECORDS')
collection = list()
for item in sliced_lst:
    z= 0
    while item[z] not in vd.keys():
        z += 1
    x = len(vd[item[z]])
    LCB = list()
    for i in range(len(item)):
        if 'empty' not in item[i]:
            z = vd[item[i]]
            z.id = str(get_first_number(item[i]))
            LCB.append(z)
        else:
            LCB.append(SeqRecord(Seq('-'*x), id=str(get_first_number(item[i]))))
    LCB = SeqIO.to_dict(LCB)
    LCB = [LCB[str(k)] for k in range(1,ngenomes+1)]
    
    # sort list of recs according to IDs
   # LCB =[r for r in sorted(LCB.values(), key=operator.attrgetter('id'))]
    
    collection.append(LCB)

def get_id(rec):
    idy = re.findall('\/([a-z,A-Z,_,0-9]*.[0-9].fasta)', rec.description)[0]
    return idy

# MERGE
print('6. MERGE')
metaaln = list()
for i in range(ngenomes):
    metarec = SeqRecord(Seq(''))
    for col in collection:
        metarec += col[i]
    metarec.id = str(i+1)
    metaaln.append(metarec)

# COLLECT PROPER IDs
print('7. COLLECT PROPER IDs')
name_dict = dict()
for key in vd.keys():
    name_dict[get_first_number(key)] = get_id(vd[key])

# ADD PROPER IDs to METAALN
print('8. ADD PROPER IDs to METAALN')
for rec in metaaln:
    rec.id = name_dict[int(rec.id)]

SeqIO.write(metaaln, outname, 'fasta')
print('FINISHED!')





