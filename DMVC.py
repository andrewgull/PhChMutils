#!/usr/bin/env python3

# a big pipeline to
# 1) take a bunch of sra reads to dump;
# 2) then to map them;
# 3) remove reads
# 4) call variants
# 5) remove BAMs

import sys
import os
import glob
import dask
from multiprocessing.pool import ThreadPool
import time


help_message = "DMVC (Download-Map-Call)- a pipeline to get VCFs starting from fastq-dumping.\n" \
               "It removes reads and BAMs to keep disk space.\n" \
               "Required arguments:\n" \
               "1 - list of samples' common names; one per line\n" \
               "2 - number of threads\n" \
               "3 - reference genome\n" \
               "4 - reference genome index name\n" \
               "5 - min alternate fraction for FreeBayes (decimal format)\n" \
               "6 - look for indels? (1/yes - yes; 0/no - no)\n" \
               "7 - fastq-dump? (1/yes - yes; 0/no - no)\n" \
               "8 - read type (p/paired - for paired reads, u/U - for unpaired)\n" \
               "9 - remove FASTQ and BAM files? (y/yes/clean/remove, n/no)\n" \
               "Example command:\n" \
               "DMVC.py samples.txt 6 H37Rv.fna H37Rv.idx 0.9 yes yes p no"

# check args
if len(sys.argv) < 7:
    print(help_message)
    sys.exit()

# args
samples_file = sys.argv[1]  # samples list
threads = int(sys.argv[2])  # number of threads to use
reference = sys.argv[3]  # reference genome name
ref_index = sys.argv[4]  # bowtie2 index base name
alt_fraction = sys.argv[5]  # alt fraction decimal
indels = sys.argv[6]  # call indels or not
dump = sys.argv[7]  # dump or not
read_type = sys.argv[8]  # paired or unpaired
clean = sys.argv[9]  # remove fastq and bam or not


def dumper(read_file):
    cmd = 'fastq-dump --split-files --clip --gzip %s' % read_file
    os.system(cmd)


def bowtie(read_file, cpus, index, paired):

    if paired == "p" or paired == "paired":
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


def freebayes(sample, frac, ref, ind):
    """
    function to run FreeBayes
    :param sample: a name of a BAM file (without extension)
    :param frac: min alternate fraction (in decimal format)
    :param ref: reference genome
    :param ind: report indels (boolean)
    :return: 0
    """
    if ind:
        print("sample %s is being processed... indels mode" % sample)
        cmd = "freebayes -f %s -p 1 --min-base-quality 20 --min-alternate-fraction %s --max-complex-gap 0 " \
              "--haplotype-length -1 --min-coverage 10 %s.bam > %s.fb.vcf" % (ref, frac, sample,  sample)
        os.system(cmd)
    else:
        print("sample %s is being processed... no indels mode" % sample)
        cmd = "freebayes -f %s -p 1 --min-base-quality 20 --min-alternate-fraction %s --max-complex-gap 0 " \
              "--haplotype-length -1 --no-indels --min-coverage 10 %s.bam > %s.fb.vcf" % (ref, frac, sample, sample)
        os.system(cmd)
    return 0


def check_file(ext):
    if len(glob.glob(ext)) == 0:
        print("Files %s not found. Exit!")
        sys.exit(1)
    else:
        pass


# read file names
samples = [line.rstrip() for line in open(samples_file).readlines()]
start_time = time.time()

# parallel fastq-dump
if dump == 1 or dump == "yes":
    print("1) start dumping from SRA-archive...")
    dask.set_options(pool=ThreadPool(threads))
    dumper_lst = list()
    for f in samples:
        dumper_lst.append(dask.delayed(dumper)(read_file=f))

    dumping_process = dask.delayed(dumper_lst)
    dumping_process.compute()

    # remove cache
    print('removing SRA cache...')
    cache_path = '/data5/bio/MolGenMicro/ncbi_cache/sra/'
    os.system('rm %s*.sra' % cache_path)
else:
    print("fastq-dump skipped....")

check_file("*.fastq.gz")

# start mapping
print("2) Mapping...")
if not os.path.isfile(ref_index + ".1.bt2"):
    os.system("bowtie2-build % %" % (reference, ref_index))
else:
    print("Genome index found...")

print("Found %i samples..." % len(samples))
for S in samples:
    bowtie(S, threads, ref_index, paired=read_type)
    samtools(S, threads)

check_file("*.bam")

# remove reads
if clean == "y" or clean == "yes" or clean == "remove":
    print("Removing FASTQ files...")
    os.system("rm *fastq.gz")

# Variant calling
print("3) Variant calling with FreeBayes...")

# for older versions of dask
dask.set_options(pool=ThreadPool(threads))
# for newer versions
# dask.config.set(pool=ThreadPool(threads))
total_lst = []
if indels == "1" or indels == "yes" or indels == "y":
    for s in samples:
        total_lst.append(dask.delayed(freebayes)(sample=s, frac=alt_fraction, ref=reference, ind=True))
elif indels == '0' or indels == 'no' or indels == 'n':
    for s in samples:
        total_lst.append(dask.delayed(freebayes)(sample=s, frac=alt_fraction, ref=reference, ind=False))

total = dask.delayed(total_lst)
total.compute()
t = time.time() - start_time

print("Finished in %s sec / %s m / %s h" % (str(t), str(t/60), str(t/3600)))

# remove BAMs
if clean == "y" or clean == "yes" or clean == "remove":
    print("removing BAM files...")
    os.system('rm *.bam')
    os.system('rm *.bam.bai')

# FINISH
print("Finished!")
