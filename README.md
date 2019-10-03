# PhChMutils
Various scripts for my projects in Physical-Chemical Medicine Institution.

Dependencies:

   - Python:
   
     1. dask
     2. pandas
     3. numpy
     4. biopython
     5. rpy2
     
   - R:
     1. VariantAnnotation
     2. GenomicFeatures
     3. data.table
     4. dplyr
     5. tidyr
     6. GenomicRanges
     7. IRanges
     8. stringr
     
   - standalone tools:
     1. samtools
     2. bowtie2
     3. FreeBayes
     4. bgzip
     5. tabix
     6. blastn
     
`ABresist.py`, `makeAnnotation.py` require all the R functions to be placed in `/home/$USER/bin`

Description of the scripts:

- `FB_parallel_dask.py` - to run FreeBayes in parallel using *dask*.
                          Run it without any args to see the help message.

- `makeAnnotation.R` - runs VariantAnnotation R package on vcf files and produces pretty table with variants and amino acid changes in each sample.
                       Works with every VCF file in a directory. Note: It strongly **dislikes** empty VCF files! Also, see *makeAnnotation.py* wrapper script.

- `makeAnnotation.py` - it indexes VCF files in a CWD and starts *VariantAnnotation.R* function then. Run it without any args to see the help message.

- `ABresist.py` - it annotates AA-changes (using *makeAnnotation.R* again) associated with AB-resistance and creates a summary tables (see below).
                  It requires additional file *resistance_SNPs_withoutrpoAC.xlsx*.
                  Run it without arguments to see the help message.
                  Table examples:
                  
                  1.Table of AB-resistance signatures in samples
 
                  | sample | INH    | RIF  | PAS  |  FQ  |
                  |------- |--------|------|------|------|
                  | ERR000 |    +   |   +  |   +  |      |
                  | ERR001 |        |      |   +  |   +  |

                  2.Table of AA changes in rpoA, rpoB and rpoC genes

                  | sample | rpoA   | rpoB  | rpoC |
                  | -------|--------|-------|------|
                  | ERR000 | A5000C | G300C |      |
                  |        | F300Y  |       |      |

- `DMVC.py` - under development!

- `MUMmerSNPs2VCF.py` - reformat *snps files produced by nucmer show-snps to VCF.

- `PhyloTyper.py` - a script for assigning each sample to phylogentic lieage using specific SNPs.
                    Also, it removes samples that have contradicting lineage signatures and writes them to a text file.
                    Run it without args to see full help message.

- `makeFASTA.py` - a script that makes pseudo-alignment file from all the VCF files in a directory.
                   It filters out PE-PPE and AB-resistance genes and invariant sites as well.
                   Run it without args to see full help message.

- `SPAtyper.py` - a script for SPA-typing (requires blast).
                  Run it without args to see full help message.
                  
- `mapping.py` - runs bowtie2 for to map your reads and samtools to process SAM and BAM files.
                 Run it without args to see full help message.

- `bgzip_and_tabix.sh` - bgzip and tabix all the VCF in a directory.

- `sff2fq.R` - very short script to convert SFF files to FASTQ (implies a certain level of your laziness).

- `collineariser.py` - a script to make meta alignment from MAUVE output files.
                       Useful for making split-graphs of bacteriophages' genomes.
                       Run it without args to see full help message.