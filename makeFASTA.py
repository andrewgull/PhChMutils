#!/usr/bin/env python3

# make fasta file from SNPs
# input: all vcf files
# output: cleaned fasta

import sys
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# from rpy2.robjects.packages import STAP
import pandas as pd
import numpy as np
import glob


def is_invariant(string):
    a = string.lower().count("a")
    t = string.lower().count("t")
    g = string.lower().count("g")
    c = string.lower().count("c")
    l = [a, t, g, c]
    if l.count(0) == 3:
        return True
    else:
        return False


def remove_invariant_sites(input_align):
    """
    removes invariant sites from an alignment
    :param input_align: Seq object = fasta of SNPs
    :return: cleaned fasta
    """
    inv = [SeqRecord(Seq('', s.seq.alphabet), id=s.id, description=s.description) for s in input_align]
    inv = MultipleSeqAlignment(inv)
    print("input alignment has %i columns" % input_align.get_alignment_length())

    for i in range(input_align.get_alignment_length()):
        if not is_invariant(input_align[:, i]):
            # add invariant column to alignment alig[:,i:i+1]
            inv = inv + input_align[:, i:i + 1]

    print("edited alignment has %i columns" % inv.get_alignment_length())
    return inv


help_message = "Script for making fasta from vcf and NJ tree from the fasta. Works on each VCF in CWD.\n" \
               "All joined vcf (result of the longest process) will be written to a TSV file\n" \
               "It uses R for tree building [turned off]\n" \
               "Five arguments:\n" \
               "1 - xlsx table with 'bad genes' and 'resistance' like ~/Documents/myco/bad_genes.xlsx\n" \
               "2 - xlsx with another bad genes, like ~/Documents/myco/bad_pos_genomes.xlsx\n" \
               "3 - gff file with reference annotation, like 'H37Rv.gff'\n" \
               "4 - name for the output fasta (a tree will be named by adding '_NJ_pdist.tre' " \
               "suffix to the fasta name)\n" \
               "5 - reference genome fasta, like H37Rv.fna\n"

if len(sys.argv) < 6:
    print(help_message)
    sys.exit()

bad_xlsx = sys.argv[1]
another_bad = sys.argv[2]
gff_name = sys.argv[3]
fasta_out = sys.argv[4]
ref_genome = sys.argv[5]

try:
    x = open(ref_genome)
except FileNotFoundError:
    print("ref genome not found!")
    sys.exit()

resistance = pd.read_excel(bad_xlsx, sheet_name="Resistance")
resistance = resistance.drop(['Gene','AA_exchange','KvarQ','Microbe_Predictor_TB','PhyResSE','TBProfiler'], axis=1)
bad_genes = pd.read_excel(bad_xlsx, sheet_name="bad_genes")
another_bads = pd.read_excel(another_bad)
another_bads = another_bads.drop(["genomes"], axis=1)
# to infer coords of bad genes you need a gff file
gff = pd.read_csv(gff_name, delimiter="\t", header=None, skiprows=7, skipfooter=1, engine="python")

# read each vcf
vcfs = glob.glob("*.vcf")
print("found %i vcf files" % len(vcfs))
print("start from %s" % vcfs[0])
concat = pd.read_csv(vcfs[0], delimiter="\t", comment="#", names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'unknown'])

concat = concat.drop(['CHROM', 'REF', 'ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'unknown'], axis=1)
# remove ALT with length more than 1 (indels)
concat = concat[concat['ALT'].map(len) < 2]
# rename ALT by sample name
concat = concat.rename(index=str, columns={"ALT": vcfs[0][:-7]})
print("concatenating vcf files...")
n = 1
for i in range(1, len(vcfs)):
    print("%i files added..." % n)
    try:
        df = pd.read_csv(vcfs[i], delimiter="\t", comment="#", names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                                                                      'FILTER', 'INFO', 'FORMAT', 'unknown'])

        # drop useless columns
        df = df.drop(['CHROM', 'REF','ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'unknown'], axis=1)
        # drop rows with ALT len more than 1
        df = df[df['ALT'].map(len) < 2]
        # rename ALT column
        df = df.rename(index=str, columns={"ALT":vcfs[i][:-7]})
        concat = pd.merge(concat, df, on="POS", how="outer")
        n += 1
    except pd.errors.EmptyDataError:
        print("empty file %s" % vcfs[i])

# sort concatenated tables by POS
concat = concat.sort_values(by="POS")
print("filtering AB resistance SNPs...")
# append resistance using set_index + join
concat = concat.set_index("POS").join(resistance.set_index("Genome_position"))
# after this operation column 'POS' turns to index, reverse it:
concat["POS"] = concat.index
# leave only NaNs in 'Antibiotic column', e.g. remove antibiotic associated positions
concat = concat[pd.isnull(concat["Antibiotic"])]

# remove another_bads using the same method: set_index + join
concat = concat.set_index("POS").join(another_bads.set_index("positions"))
concat["POS"] = concat.index
concat = concat[pd.isnull(concat["effect"])]

# clean concat from bad genes etc
# start from bad genes -  clean by gff
PPE = bad_genes["Locus tag"]
# leave gff with only PPE genes
gff = gff[gff[8].map(lambda x: any(gene in x for gene in PPE))]
PPE_coords = gff[[3,4]]
# remove PPE genes from concat
coords_tuples = [tuple(x) for x in PPE_coords.values]
print("filtering PE-PPE...")
for coord_pair in coords_tuples:
    concat = concat[concat["POS"].map(lambda x: x not in range(coord_pair[0], coord_pair[1]+1))]

# remove useless columns 'POS', 'Antibiotic' and 'effect';
concat = concat.drop(['Antibiotic', 'effect'], axis=1)

# make reference column with bases from ref corresponding positions
ref = SeqIO.read(open(ref_genome), "fasta")
ref_column = [str(ref.seq[p-1]) for p in concat["POS"]]
concat["REF"] = ref_column

# substitute NaNs in each columns (like if_else in R)
for name in list(concat):
    try:
        concat[name] = np.where(pd.isnull(concat[name]), concat['REF'], concat[name])
    except ValueError:
        print("value error - %s" % name)

# drop REF column
concat = concat.drop(["REF"], axis=1)
# write resulted concat table to file
concat.to_csv("concatenated_"+str(len(vcfs))+"vcf.tsv", sep="\t", index=False)
# collect columns as character strings
records = []
for col in list(concat)[:len(list(concat))-1]:
    snps = list(concat[col])
    snps = Seq("".join(snps))
    rec = SeqRecord(snps, id=col, description="")
    records.append(rec)


snps_inv = remove_invariant_sites(MultipleSeqAlignment(records))
AlignIO.write(snps_inv, fasta_out, "fasta")
print("Invariant SNPs are ready...")

# while Bio.Phylo.TreeConstruction is not available use R to build a NJ-tree
buildNJ_R = """
buildNJ <- function(in.fa.file){
  library(ape); library(phangorn)
  seqs <- read.dna(in.fa.file, format="fasta")
  seqs.phyDat <- phyDat(seqs, type="DNA", levels=NULL)
  dist.mat.pdist <- dist.p(seqs.phyDat)
  NJ_tree_pdist <- NJ(dist.mat.pdist)
  write.tree(NJ_tree_pdist, file=paste(in.fa.file, "NJ_pdist.tre", sep="_"))
}
"""

# buildTree = STAP(buildNJ_R, "buildTree")
# buildTree.buildNJ(fasta_out)
print("Done! NJ building was turned off")
