#!/usr/bin/env python3
# script for running FreeBayes in parallel

from multiprocessing.pool import ThreadPool
import os
import sys
import dask
import time

start_time = time.time()
help_message = "FreeBayes paralleled by dask, v. 1.03. Arguments:\n" \
               "1 - list of samples' common names; one per line\n" \
               "2 - number of threads\n" \
               "3 - min alternate fraction (decimal format)\n" \
               "4 - reference genome\n" \
               "5 - indels (1, yes, y - report indels; 0, no, n - do not report indels)"
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


# FB function
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


# run in parallel
# for older versions of dask
dask.set_options(pool=ThreadPool(threads))
# for newer versions
# dask.config.set(pool=ThreadPool(threads))
total_lst = []
if indels == "1" or indels == "yes" or indels == "y":
    for s in cmn_names:
        total_lst.append(dask.delayed(freebayes)(sample=s, frac=alt_fraction, ref=reference, ind=True))
elif indels == '0' or indels == 'no' or indels == 'n':
    for s in cmn_names:
        total_lst.append(dask.delayed(freebayes)(sample=s, frac=alt_fraction, ref=reference, ind=False))
else:
    print("'indels' argument was not recognised."
          "Use either 'yes', 'y' or '1' to report indels or 'no', 'n' or '0' to not report them")
    sys.exit()

total = dask.delayed(total_lst)
total.compute()
t = time.time() - start_time

print("Finished in %s sec / %s m / %s h" % (str(t), str(t/60), str(t/3600)))
