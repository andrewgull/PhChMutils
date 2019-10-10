#!/usr/bin/env python3

# wrapper script to run kallisto and quant and organize output files

import os
import sys


def kallisto(r1, r2, threads=10, boot=100):
    os.system('kallisto quant -i H37Rv.cds.idx -o kallisto_dir -b %s -t %s --genomebam -g H37Rv.gtf %s %s' % (boot, threads, r1, r2))


if len(sys.argv) < 5:
    print("wrapper for 'kallisto quant'.\nRef genome GTF file must be in CWD.\n"
          "Args:\n"
          "1 - reads file (one per line)\n"
          "2 - threads\n"
          "3 - bootstrap\n"
          "4 - number of characters to cut from the end of a filename to get sample name\n"
          "(e.g. 11 for 'SAMPLE_1_fastq.gz'")
    sys.exit()

readsfile = sys.argv[1]
thrds = sys.argv[2]
bt = sys.argv[3]
cut = int(sys.argv[4])

reads = [line.rstrip() for line in open(readsfile)]
# ['MT_1947_1_R1.fq.gz', 'MT_1947_1_R2.fq.gz']
for i in range(0, len(reads), 2):
    print('processing %s...' % reads[i][:-cut])
    kallisto(reads[i], reads[i+1], threads=thrds, boot=bt)
    print('renaming')
    os.system('mv kallisto_dir/abundance.tsv kallisto_dir/abundance_%s.tsv' % reads[i][:-cut])
    os.system('mv kallisto_dir/pseudoalignments.bam kallisto_dir/pseudoalignments_%s.bam' % reads[i][:-cut])
    os.system('mv kallisto_dir/pseudoalignments.bam.bai kallisto_dir/pseudoalignments_%s.bam.bai' % reads[i][:-cut])

os.chdir("./kallisto_dir")
# extract the 4th column from abundance files
bash_extract = """for F in abundance_*.tsv; do cut -f4 -d$'\t' $F | tail -n +2 > $F".tpm"; echo -e ${F:10:10} | cat - $F".tpm" > $F".head.tpm"; cut -f1 -d$'\t' $F > feature_column; paste feature_column *.head.tpm > TPM_ALL.tsv; done"""

os.system(bash_extract)

print('done')
