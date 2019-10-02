#!/usr/bin/env python3

from rpy2.robjects.packages import STAP
import glob
import os
import getpass
import sys
import warnings
from rpy2.rinterface import RRuntimeWarning
warnings.filterwarnings("ignore", category=RRuntimeWarning)


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


user_name = getpass.getuser()

# makeAnnotation.R is called by another function because it's easier to make changes
# in a single copy of makeAnnotation.R stored in /home/$USER/bin
runAnnotation = """
runAnnotation <- function(ref_genome, ref_gff, NC, species_name, with_ranges=FALSE, in_ranges=NULL, cpus=6){
  source("/home/%s/bin/makeAnnotation.R")
  annots <- makeAnnotation(ref.genome.name=ref_genome, ref.annotation=ref_gff, nc=NC, species=species_name, 
  withRanges=with_ranges, inRanges=in_ranges, writeXLSX=FALSE, outTable="SNP_annotation_table", substPart=".vcf.gz", 
  pat=".vcf", threads=cpus)
}
""" % user_name

R_Annot = STAP(runAnnotation, "R_Annot")

help_message = "Script for annotating AA-changes in VCF files\n" \
               "VCF files are taken from the CWD\n" \
               "makeAnnotation.R function should be in ~/bin\n" \
               "Required arguments:\n" \
               "1 - reference genome accession (NC_000962)\n" \
               "2 - reference genome fasta (H37Rv.fna)\n" \
               "3 - reference genome annotation (H37Rv.gff)\n" \
               "4 - species name ('Mycobacterium tuberculosis')\n" \
               "5 - number of threads\n"

if len(sys.argv) < 6:
    print(help_message)
    sys.exit()

fun_path = "/home/%s/bin/makeAnnotation.R" % user_name
while not os.path.isfile(fun_path):
    print("File %s not found!" % fun_path)
    fun_path = input("Enter full path to the 'makeAnnotation.R' function > ")

vcf_files = sorted(glob.glob("*.vcf"))
genome_nc = sys.argv[1]
genome_fa = sys.argv[2]
genome_gff = sys.argv[3]
species = sys.argv[4]
threads = sys.argv[5]

# bgzip and tabix VCF files
if len(glob.glob("*.tbi")) != len(vcf_files):
    print("vcf bgzipping and tabixing...")
    bgzip_and_tabix(vcf_files)
else:
    print("Found indexed VCF files, skip indexing...")

# check the genome index
if not os.path.isfile(genome_fa + ".fai"):
    print("Genome index not found...\nRun samtools faidx...")
    os.system("samtools faidx %s" % genome_fa)
else:
    print("Found genome index, skip indexing...")

# run makeAnnotation.R
print("Annotation...")
R_Annot.runAnnotation(genome_fa, genome_gff, genome_nc, species, cpus=threads)
