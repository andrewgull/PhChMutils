#!/usr/bin/env python3

# wrapper script to run kallisto quant and organize output files

import os
import sys
import glob
from Bio import SeqIO
import re
import pandas as pd


def kallisto(r1, r2, index="H37Rv.cds.idx", threads=1, boot=10):
    os.system('kallisto quant -i %s -o kallisto_dir -b %s -t %s --genomebam -g H37Rv.gtf %s %s' % (index, boot, threads, r1, r2))


if len(sys.argv) < 5:
    print("wrapper for 'kallisto quant'.\nRef genome GTF file and indexed CDS/transcripts (*.idx) must be in CWD.\n"
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

index_file = glob.glob("*.idx")[0]
gtf_file = glob.glob("*.gtf")[0]
if not os.path.isfile(index_file) or not os.path.isfile(gtf_file):
    print("Index or GTF file not found!\nTo make index run 'kallisto index -i index_name.idx transcripts_file.fasta'")
    sys.exit(1)

reads = [line.rstrip() for line in open(readsfile)]
# ['MT_1947_1_R1.fq.gz', 'MT_1947_1_R2.fq.gz']
for i in range(0, len(reads), 2):
    print('processing %s...' % reads[i][:-cut])
    kallisto(reads[i], reads[i+1], index=index_file, threads=thrds, boot=bt)
    print('renaming')
    os.system('mv kallisto_out/abundance.tsv kallisto_dir/abundance_%s.tsv' % reads[i][:-cut])
    os.system('mv kallisto_out/pseudoalignments.bam kallisto_dir/pseudoalignments_%s.bam' % reads[i][:-cut])
    os.system('mv kallisto_out/pseudoalignments.bam.bai kallisto_dir/pseudoalignments_%s.bam.bai' % reads[i][:-cut])

os.chdir("./kallisto_out")
# extract the 4th column from abundance files
# original bash version
# for F in abundance_*.tsv;
#   do cut -f4 -d$'\t' $F | tail -n +2 > $F".tpm";
#   echo -e ${F:10:10} | cat - $F".tpm" > $F".head.tpm";
#   cut -f1 -d$'\t' $F > feature_column;
#   paste feature_column *.head.tpm > COUNTS_ALL.tsv;
# done

# python version
init_tab = pd.read_csv(glob.glob("abundance_*.tsv")[0], delimiter="\t")[['target_id']]
init_tab = init_tab.set_index('target_id')

for f in glob.glob("abundance_*.tsv"):
    table = pd.read_csv(f, delimiter="\t", index_col=None)
    table = table[['target_id', 'est_counts']]
    table.columns = ['target_id', f[10:-4]]
    table = table.set_index('target_id')
    init_tab = table.join(init_tab)

# add column with loci tags and gene names
cds_filename = "../" + index_file[:-3] + "fasta"
cds_headers = [rec.description for rec in SeqIO.parse(cds_filename, "fasta")]
ids = [rec.id for rec in SeqIO.parse(cds_filename, "fasta")]
loci_tags = list()
gene_names = list()
for rec in cds_headers:
    try:
        gene = re.findall(r'\[gene=(.*)\] \[locus_tag=.*\]', rec)[0]
    except IndexError:
        gene = ''
    try:
        tag = re.findall(r'\[locus_tag=(.*)\] \[db_xref=.*\]', rec)[0]
    except IndexError:
        tag = ''
    gene_names.append(gene)
    loci_tags.append(tag)

names_and_tags = pd.DataFrame({'target_id': ids, 'gene': gene_names, 'locus_tag': loci_tags})

full_table = names_and_tags.set_index('target_id').join(init_tab)

# write to a file
full_table.to_csv('counts_all.tsv', sep='\t', na_rep="")

print('done')
