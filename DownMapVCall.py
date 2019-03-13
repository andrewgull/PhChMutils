#! /srv/common/bin/python3
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


# read a list of files to dump
help_message = "DownMapVCall - a pipeline to get VCFs starting from SRA-dumping. " \
               "Removes reads and BAMs to keep disk space. " \
               "Arguments:\n1 - list of samples' common names; one per line\n" \
               "2 - number of threads\n" \
               "3 - reference genome\n" \
               "4 - min alternate fraction for FreeBayes (decimal format)\n" \
               "5 - look for indels? (1 - yes; 0 - no)\n" \
               "6 - fastq-dump? (1 -yes; 0 - no)"

# check args
if len(sys.argv) < 5:
    print(help_message)
    sys.exit()

# check file with names
try:
    filename = sys.argv[1]
except FileNotFoundError:
    print("Error! File of samples not found")
    sys.exit()

# args
threads = int(sys.argv[2])
reference = sys.argv[3]
alt_fraction = sys.argv[4]
indels = sys.argv[5]
dump = sys.argv[6]

# read file names
with open(filename) as f:
    cmn_names = f.readlines()
cmn_names = [x.strip() for x in cmn_names]


######### FUNCTIONS ######################
# dump function
def dumper(read_file):
    cmd = 'fastq-dump --split-files --clip --gzip %s' % read_file
    os.system(cmd)


def bowtie(read_file, cpus, index):
    r1 = read_file + '_1.fastq.gz'
    r2 = read_file + '_2.fastq.gz'
    bowtie_cmd = 'bowtie2 -q -t --local -x %s -p %s -1 %s -2 %s  -S  %s.sam 2> %s.stats' % (index, cpus, r1, r2,
                                                                                            read_file, read_file)
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


# FB fun (with indels)
def freebayes(sample, frac, ref):
    print("sample %s is being processed... indels mode" % sample)
    cmd = "freebayes -f %s -p 1 --min-base-quality 20 --min-alternate-fraction %s --max-complex-gap 0 " \
          "--haplotype-length -1 --min-coverage 10 %s.bam > %s.fb.vcf" % (ref, frac, sample, sample)
    os.system(cmd)
    return 0


# FB fun (no indels)
def freebayes_noind(sample, frac, ref):
    print("sample %s is being processed... no indels mode" % sample)
    cmd = "freebayes -f %s -p 1 --min-base-quality 20 --min-alternate-fraction %s --max-complex-gap 0 " \
          "--haplotype-length -1 --no-indels --min-coverage 10 %s.bam > %s.fb.vcf" % (ref, frac, sample, sample)
    os.system(cmd)
    return 0


########## THE PIPELINE ################

# parallel fastq-dump
if dump == 1 or dump == "yes":
    print("start dumping from SRA-archive...")
    dask.set_options(pool=ThreadPool(threads))
    dumper_lst = list()
    for f in cmn_names:
        dumper_lst.append(dask.delayed(dumper)(read_file=f))

    dumping_process = dask.delayed(dumper_lst)
    dumping_process.compute()

    # remove cache
    print('removing SRA cache...')
    cache_path = '/data5/bio/MolGenMicro/ncbi_cache/sra/'
    os.system('rm %s*.sra' % cache_path)
else:
    print("No dumping was selected. Go to mapping.")

# start mapping
# build index
if ".fna" in reference:
    index_name = reference[:-4]
elif ".fasta" in reference:
    index_name = reference[:-6]
elif ".fa" in reference:
    index_name = reference[:-3]
else:
    print("Reference genome file extension is not recognized. It must be '.fa', '.fna' or '.fasta'.\nQuit.")
    sys.exit()

if len(glob.glob("*.bt2")) == 6:
    print("It seems there are index files in the CWD. Skip genome indexing.")
else:
    print("No genome index found. Indexing...")
    os.system('bowtie2-build %s %s' % (reference, index_name))

print("mapping...")
print("there are %i pairs of files" % len(cmn_names))
n = 0
for f in cmn_names:
    print("working on %s..." % f)
    bowtie(f, cpus=threads, index=index_name)
    samtools(f, cpus=threads)
    n += 1
    if n % 250 == 0:
        print('\n%i pairs are mapped\n' % n)

# remove reads
print("removing FASTQ files...")
os.system("rm *fastq.gz")

# Variant calling
print("Variant calling with FreeBayes in parallel...")

# run in parallel
dask.set_options(pool=ThreadPool(threads))
total_lst = []
if indels == "1":
    # with indels
    for s in cmn_names:
        total_lst.append(dask.delayed(freebayes)(sample=s, frac=alt_fraction, ref=reference))
else:
    for s in cmn_names:
        total_lst.append(dask.delayed(freebayes_noind)(sample=s, frac=alt_fraction, ref=reference))

total = dask.delayed(total_lst)
total.compute()

# remove BAMs
print("removing BAM files...")
os.system('rm *.bam')
os.system('rm *.bam.bai')

# FINISH
print("dumping-mapping-calling done!")
