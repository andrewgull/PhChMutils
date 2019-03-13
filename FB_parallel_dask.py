#!/usr/bin/python3
# script for running FreeBayes in parallel

from multiprocessing.pool import ThreadPool
import os
import sys
import dask
import time

start_time = time.time()
help_message = "FreeBayes paralleled by dask. Arguments:\n1 - list of samples' common names; one per line\n2 - number " \
               "of threads\n3 - min alternate fraction (decimal format)" \
               "\n4 - reference genome\n5 - indels (1 - yes; 0 - no)"
if len(sys.argv) < 4:
    print(help_message)
    sys.exit()

try:
    filename = sys.argv[1]
except FileNotFoundError:
    print("Error! File of samples not found")
    sys.exit()

threads = int(sys.argv[2])
alt_fraction = sys.argv[3]
reference = sys.argv[4]
indels = sys.argv[5]

# read file with common names
with open(filename) as f:
    cmn_names = f.readlines()
cmn_names = [x.strip() for x in cmn_names]


# FB fun (with indels)
def freebayes(sample, frac, ref):
    print("sample %s is being processed... indels mode" % sample)
    cmd = "freebayes -f %s -p 1 --min-base-quality 20 --min-alternate-fraction %s --max-complex-gap 0 " \
          "--haplotype-length -1 --min-coverage 10 %s.bam > %s.fb.vcf" % (ref, frac, sample,  sample)
    os.system(cmd)
    return 0


# FB fun (no indels)
def freebayes_noind(sample, frac, ref):
    print("sample %s is being processed... no indels mode" % sample)
    cmd = "freebayes -f %s -p 1 --min-base-quality 20 --min-alternate-fraction %s --max-complex-gap 0 " \
          "--haplotype-length -1 --no-indels --min-coverage 10 %s.bam > %s.fb.vcf" % (ref, frac, sample,  sample)
    os.system(cmd)
    return 0


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

print("Finished in %s sec" % (time.time() - start_time))
