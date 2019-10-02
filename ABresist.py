#!/usr/bin/env python3

# script for AB resistance
# it should transform vcf with 'transform_vcf.sh'
# and make all the steps to run R function 'makeAnnotation.R' (including IRanges object)
# everything works fine, but intergenic regions for all genome are in the resulting table
# because they're taken from vcf files, not from GRanges regions


from rpy2.robjects.packages import STAP
import pandas as pd
import glob
import os
import re
import sys
from shutil import copy
import csv
from collections import defaultdict


def bgzip_and_tabix(files):
    """
    fun to bgzip and tabix transformed VCFs
    :param files: vcf files to bgzip an tabix
    :return: no
    """
    for f in files:
        os.system("bgzip %s" % f)
        os.system("tabix -p vcf %s" % (f+".gz"))
    print("bgzipping and tabixing is finished")


def get_column(csv_tab, column):
    """
    function to extract column from csv tab (list)
    :param csv_tab: table as list of rows
    :param column: column to be extracted
    :return: column as a list
    """
    idx = csv_tab[0].index(column)
    col = [row[idx] for row in csv_tab]
    return col


user_name = os.getlogin()

# here goes a function to make GRanges object of AB-resistance genes
resist2GRanges = """
resist2GRanges <- function(resdf){
  # a function to make IRanges object of AB-resistance genes 
  # first three column names: "Genome position", "Gene", "AA exchange"
  
  # inRanges - GRanges object of regions of interest
  # example: GRanges(seqnames = "NC_000962", ranges = IRanges(start=c(3877464, 759807, 763370), end=c(3878507, 763325, 767320), names=c("rpoA", "rpoB", "rpoC")))
  
  resdf <- resdf[,c(1:3)]
  genes <- unique(resdf[,2])
  res.list <- split(resdf, resdf$Gene)
  gene.starts <- sapply(res.list, function(x){min(x[,1])})
  gene.ends <- sapply(res.list, function(x){max(x[,1])})
  range_obj <- inRanges <- GRanges(seqnames="NC_000962", ranges=IRanges(start=gene.starts, end=gene.ends, names=sort(genes$Gene)))
  return(range_obj)
}
"""

# a function uniting all other R functions need for annotation
fullAnnotation = """
fullAnnotationInGRanges <- function(resistance_table){
  source("/home/%s/bin/makeAnnotation.R")
  source("/home/%s/bin/resist2GRanges.R")
  resist_ranges <- resist2GRanges(resistance_table)
  annots <- makeAnnotation(withRanges=TRUE, inRanges=resist_ranges, outTable="annotation_rpoABC")
  
}
""" % (user_name, user_name)

R_Annot = STAP(fullAnnotation, "R_Annot")

help_message = "Script for annotating AA-changes in AB-resistance regions, v 1.2\n" \
               "VCF files are taken from the CWD\n" \
               "Resistance is taken from resistance_SNPs_withoutrpoAC.csv table\n" \
               "Three R functions should be in ~/bin\n" \
               "Four arguments have to be passed:\n" \
               "1 - reference genome accession (NC_000962)\n" \
               "2 - reference genome fasta (H37Rv.fna)\n" \
               "3 - reference genome annotation (H37Rv.gff)"

if len(sys.argv) < 4:
    print(help_message)
    sys.exit()

# check how many vcf files have SNPs in resistance table
# read table and files and check whether they have AB mutations
r_tab = "/data5/bio/MolGenMicro/mycobacterium/resistance_SNPs_withoutrpoAC.csv"

if not os.path.isfile(r_tab):
    print("Resistance file %s not found" % r_tab)
    r_tab = input("Enter full path to the file > ")

vcf_files = sorted(glob.glob("*.vcf"))
genome_nc = sys.argv[1]
genome_fa = sys.argv[2]
genome_gff = sys.argv[3]

# read resistance table and collect genomic positions from there
# TAKE CARE of the footer
res_tab = pd.read_csv(r_tab, nrows=1370, delimiter='\t')
mut_pos = res_tab["Pos"]

r_functions = ["/home/%s/bin/makeAnnotation.R" % user_name, "/home/%s/bin/resist2GRanges.R" % user_name,
               "/home/%s/bin/fullAnnotationInGRanges.R" % user_name]

for fun in r_functions:
    if not os.path.isfile(fun):
        print("file %s not found" % fun)
        sys.exit(1)

# bgzip and tabix VCF files
if len(glob.glob("*.tbi")) != len(vcf_files):
    print("vcf bgzipping and tabixing...")
    bgzip_and_tabix(vcf_files)
else:
    print("Found indexed VCF files, skip indexing...")

# index the genome
if len(glob.glob("*.fai")) == 0:
    os.system("samtools faidx %s" % genome_fa)

print("Annotation...")
R_Annot.fullAnnotationInGRanges(r_tab)

# Table with resistance
#
# sample | INH    | RIF  | PAS  |  FQ  |
# -------|--------|------|------|------|
# ERR000 |    +   |   +  |   +  |      |
#--------|--------|------|------|------|
# ERR001 |        |      |   +  |   +  |
#--------|--------|------|------|------|
#

print("making resistance table...")

# R function to join annotation and AB resistance tables
join_tables = """
join_tables <- function(path1, path2){
  # function to join tables: path1 = path to AB resist - "~/Documents/myco/compensatory/resistance_SNPs_withoutrpoAC.csv";
  # path2 = path to annotation table - "./annotation_table.tsv"
  library(dplyr)
  resistanceSNP <- read.csv(path1, nrows=1370, sep='\t')
  resistanceSNP <- resistanceSNP[,c(1:4)]
  annotation_table <- read.delim(path2)
  # join your annotation with resistance table
  both <- left_join(annotation_table, resistanceSNP, by="Pos")
  both <- select(both, c(Gene, `AA exchange`, Antibiotic), everything())
  #n <- colnames(both)[4:16]
  both <- select(both, c(Chrom, Pos, Ref, strand, CDSLOC.start, CDSLOC.end, PROTEINLOC, GENEID, CONSEQUENCE, REFCODON, VARCODON, REFAA, VARAA), everything())
  write.csv(both, "annotation+resistance.csv", row.names=F, quote=F, na="")
}
"""

# load the function
joiner = STAP(join_tables, "joiner")

# and run it
print("joining tables...")
joiner.join_tables(r_tab, "annotation_rpoABC.tsv")

# after joined table is written you can read it and iterate

# read annotation+resistance table
# remove rows with no AB
with open("annotation+resistance.csv") as csvfile:
    annotation = csv.reader(csvfile, delimiter=",")
    annotation = [row for row in annotation]
# leave only with antibiotics
annotation = [row for row in annotation if row[15] != ""]
samples = sorted(annotation[0][16:])

# for sample in samples:
#   for row in rows of annotation:
#     if in cell stands a letter, then "+" to corr AB key in resist_dict

# get non redundant list of all ABs
antibiotic_names = set(res_tab["Antibiotic"])
# make blank dict (=table)
resist_dict = defaultdict(list)
resist_dict["samples"] = samples
for ant in antibiotic_names:
    resist_dict[ant] = [""] * len(samples)

# iterate through samples and leave "+" if AB resistance
print("started iteration to make table with resistance marks...")
for j in range(len(samples)):
    sample_column = get_column(annotation, samples[j])
    for i in range(1,len(annotation)):
        if sample_column[i] != "":
            resist_dict[annotation[i][15]][j] = "+"

# make DataFrame
resist_table = pd.DataFrame.from_dict(resist_dict)

# Table with SNPs and AA changes in rpoA, B, C
#
# sample | rpoA   | rpoB  | rpoC |
# -------|--------|-------|------|
# ERR000 | A5000C | G300C |
#        | F300Y

# read vcf file and if an SNP falls between rpoABC coords
# add it to special table
rpoA = (3877464,3878507)  # complement strand
rpoB = (759807,763325)
rpoC = (763370,767320)

# make such blank table (=dict)
rpo_dict = defaultdict(list)

for v in vcf_files:
    print(v)
    with open(v) as f:
        vcf = f.readlines()
    # get rid of the header
    vcf = [line for line in vcf if "#" not in line]
    snp_positions = [re.findall("NC_000962\t([0-9]*)\t.", line)[0] for line in vcf]
    substitutions = [re.findall("NC_000962\t([0-9]*\t.\t[A,T,G,C]*\t[A,T,G,C]*)\t", line)[0].split("\t") for line in vcf]
    substitutions = ["-".join(item[-2:]+[item[0]]) for item in substitutions]
    for p in snp_positions:
        idx = snp_positions.index(p)
        if int(p) in range(rpoA[0], rpoA[1]):
            rpo_dict["rpoA"].append(substitutions[idx])
            rpo_dict["samples"].append(v)
            rpo_dict["rpoB"].append("")
            rpo_dict["rpoC"].append("")
        elif int(p) in range(rpoB[0], rpoB[1]):
            rpo_dict["rpoB"].append(substitutions[idx])
            rpo_dict["samples"].append(v)
            rpo_dict["rpoA"].append("")
            rpo_dict["rpoC"].append("")
        elif int(p) in range(rpoC[0], rpoB[0]):
            rpo_dict["rpoC"].append(substitutions[idx])
            rpo_dict["samples"].append(v)
            rpo_dict["rpoA"].append("")
            rpo_dict["rpoB"].append("")

# Only AA changes are absent :(
rpo_table = pd.DataFrame.from_dict(rpo_dict)
rpo_table["samples"] = [s[:-4] for s in rpo_table["samples"]]

# now join two tables by sample
concat = resist_table.set_index("samples").join(rpo_table.set_index("samples"))
concat.to_csv("concat_tab.csv", index=False)
