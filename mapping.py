#!/usr/bin/env python3

import os
import sys

help_message = "Script for mapping reads with bowtie2 and processing with samtools v1.9\n" \
               "Four arguments have to be passed:\n" \
               "1 - paired or not? (p/P - for paired reads, u/U - for unpaired)\n" \
               "2 - file with samples' names, one per line\n" \
               "3 - genome index base name (use 'bowtie2-build genome index_name')\n" \
               "4 - number of threads"

if len(sys.argv) < 5:
    print(help_message)
    sys.exit()

read_type = sys.argv[1]
samples_file = sys.argv[2]
ref_index = sys.argv[3]
threads = sys.argv[4]


def bowtie(read_file, cpus, index, paired=True):

    if paired:
        r1 = read_file + '_1.fastq.gz'
        r2 = read_file + '_2.fastq.gz'
        bowtie_cmd = 'bowtie2 -q -t --local -x %s -p %s -1 %s -2 %s  -S  %s.sam 2> %s.stats' % (index, cpus, r1, r2,
                                                                                            read_file, read_file)
    else:
        r1 = read_file + '.fastq.gz'
        bowtie_cmd = 'bowtie2 -q -t --local -x %s -p %s -U %s -S  %s.sam 2> %s.stats' % (index, cpus, r1, read_file,
                                                                                         read_file)
    os.system(bowtie_cmd)


def samtools(read_file, cpus):
    view_cmd = "samtools view -bS -@ %s %s.sam > %s.bam" % (cpus, read_file, read_file)
    sort_cmd = "samtools sort -o %s.sorted.bam -O BAM -@ %s %s.bam" % (read_file, cpus, read_file)
    index_cmd = "samtools index %s.bam" % read_file

    os.system(view_cmd)
    os.system("rm %s.sam" % read_file)
    os.system(sort_cmd)
    os.system("rm %s.bam" % read_file)
    os.system("mv %s.sorted.bam %s.bam" % (read_file, read_file))
    os.system(index_cmd)


samples = [line.rstrip() for line in open(samples_file).readlines()]

if read_type == 'p' or 'P':
    print("Found %i pairs of files..." % len(samples))
else:
    print("Found %i files..." % len(samples))

for S in samples:
    if read_type == 'p' or 'P':
        bowtie(S, threads, ref_index, paired=True)
    else:
        bowtie(S, threads, ref_index, paired=False)

    samtools(S, threads)

print("Mapping finished successfully")
