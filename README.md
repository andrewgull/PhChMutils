# PhChMutils
Various scripts for my projects in Physical-Chemical Medicine Institution:

- `FB_parallel_dask.py` - runs FreeBayes in parallel using *dask*

- `makeAnnotation.R` - runs VariantAnnotation package on vcf files to produce pretty table with variants and amino acid changes

- `ABresist.py` - a script for annotating AA-changes in antibiotic resistance regions. Run it without arguments and see full help message.

- `DMVC.py` - a pipeline to:
    1. download SRA files 
    2. map reads on reference 
    3. free disk space by removing already mapped SRA 
    4. call variants with FreeBayes 
    5. free disk space by removing *bam* files

    run it without args to see full help message.

- `MUMmerSNPs2VCF.py` - reformat *snps files produced by MUMmer to VCF.

- `PhyloTyper.py` - a script for assigning each sample to phylogentic lieage using specific SNPs. Also it removes samples that got contradicting lineage signatures and writes them to text file.

- `makeFASTA.py` - a script that makes pseudo-alignment file from all the VCF files in a directory. Filters out PE-PPE and AB-resistance genes and invariant sites.

- `SPAtyper.py` - a script for SPA-typing (requires blast).

- `bgzip_and_tabix.sh`

- `snps2vcf.R`
