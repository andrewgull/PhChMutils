makeAnnotation <- function(ref.genome.name="H37Rv.fna", ref.annotation="H37Rv.gff", nc="NC_000962", species="Mycobacterium tuberculosis", withRanges=FALSE, inRanges=NULL, writeXLSX=FALSE, outTable="annotation_table", substPart=".vcf.gz", pat=".vcf", threads=8){
  # the function makes annotation table with amino acid changes in wide format for the whole genome or for a region specified;
  # for flawless work of the function the following conditions have to be true:
  # 1) all vcf files have to be bgzipped and tabixed
  # 2) reference genome have to be indexed with samtools faidx
  # 3) original VCFs have to be in the CWD along with gzipped and tabixed versions
  # rrna.genes <- vector of vectors of rRNA genes' coordinates, like c(c(1471846, 1473382), c(1473658, 1476795))
  # Deprecated: 1) all vcf files were transformed with sed -i '/##/d' and have extension 'vcf'
  

  # inRanges - GRanges object of regions of interest
  # example: GRanges(seqnames = "NC_000962", ranges = IRanges(start=c(3877464, 759807, 763370), end=c(3878507, 763325, 767320), names=c("rpoA", "rpoB", "rpoC")))
  # withRanges have to be set TRUE
  # pat <- part of a file name that have to substituted while making vcf table (sample name) and pattern to find that type of files
  # substPart <- part of a file name (gzipped vcf) that have to substituted while making prediction table
  # threads - number of threads to use in mclapply() function
  # return - a list of two tables: full annotation long and full annotation wide
  
  # source("https://bioconductor.org/biocLite.R")
  library(parallel);  library(VariantAnnotation);  library(GenomicFeatures);  library(dplyr);  library(stringr);  library(data.table);  library(tidyr)
  #library(xlsx)
  
  
  ### THE INTERNAL FUNCTIONS BLOCK ################################################################
  transform_col <- function(column){
    # transform the cell from char vector to chars
    # column arg should look like 'table$column'
    clist <- CharacterList(column)
    mult <- elementNROWS(clist) > 1L
    # for older versions of libraries or even R
    #mult <- elementLengths(clist) > 1L
    clist[mult] <- lapply(clist[mult], paste0, collapse=",")
    return(as.character(clist))
  }
  
  vcf2table <- function(files=files.vcf){
    # a function to bind all vcf files in single table
    # 
    vcf.list <- lapply(files, function(x){
      t <- fread(x, sep = "\t", header = T, select = c("#CHROM", "POS", "REF", "ALT"), skip="#CHROM")
      sample.name <- sub(pat, "", x, fixed=T)
      t$sample <- rep(sample.name, nrow(t))
      return(t)
    })
    vcf.table <- rbindlist(vcf.list)
    return(vcf.table)
  }
  
  changeInRegion <- function(fl, rng, refgenome, annotation){
    # a function to annotate amino acid changes in specific region
    tab <- TabixFile(fl)
    vcf.rng <- readVcf(tab, nc, param=rng)
    aa <- predictCoding(vcf.rng, annotation, seqSource = refgenome)
    #aa <- as.data.frame(aa)
    aa <- as.data.frame(aa, row.names = c(1:length(aa)))
    sample.name <- sub(".vcf.gz", "", fl, fixed = TRUE)
    aa$sample <- rep(sample.name, nrow(aa))
    return(aa)
  }
  
  changeAll <- function(fl, ref, anntn){
    # a function to annotate amino acid changes in whole genome sequence
    #print(fl)
    tab <- TabixFile(fl)
    vcf.rng <- readVcf(tab, nc)
    aa <- predictCoding(vcf.rng, anntn, seqSource = ref)
    #aa <- as.data.frame(aa)
    if (length(aa) == 0){
      print(fl)
      print(length(aa))
    } else {
      aa <- as.data.frame(aa, row.names = c(1:length(aa)))
      sample.name <- sub(substPart, "", fl, fixed = TRUE)
      aa$sample <- rep(sample.name, nrow(aa))
      return(aa)
    } 
  }
  
  prettyTable <- function(tab){
    # a function to make annotation table look much better
    tab$ALT <- transform_col(tab$ALT)
    tab$PROTEINLOC <- transform_col(tab$PROTEINLOC)
    #tab <- select(tab, -c(seqnames, end, paramRangeID, QUAL, FILTER, varAllele, CDSLOC.width, QUERYID, TXID, CDSID))
    tab <- tab[, -c(1, 3, 6, 9, 10, 11, 14, 16, 17, 18)]
    tab$sample <- sub(substPart, "", tab$sample, fixed=T)
    return(tab)
  }
  #### END OF THE INTERNAL FUNCTIONS BLOCK #############################################################
  
  # check if genome was indexed
  idx <- dir(pattern=".fai")
  if (length(idx)==0){
    print("Index your genome with samtools faidx first!")
    break()
  }
  
  # prepare reference annotation
  refgenome <- FaFile(ref.genome.name)
  annot <- makeTxDbFromGFF(ref.annotation, format = "auto", organism=species, circ_seqs = nc)
  
  # bgzip and tabix all vcf files in the CWD
  # Also: 'samtools faidx genome.fna'
  # get the list of them
  files.vcf.gz <- dir(pattern='.vcf.gz$')
  files.tbi <- dir(pattern=".tbi$")
  files.vcf <- dir(pattern='.vcf$')
  
  # make annotation table for a region(s)?
  # prepare the regions
  # example
  # myrange <- GRanges(seqnames = "NC_000962", ranges = IRanges(start=c(3877464, 759807, 763370), end=c(3878507, 763325, 767320), names=c("rpoA", "rpoB", "rpoC")))
  if (withRanges){
    myrange <- inRanges
    predict.list <- mclapply(files.vcf.gz, function(x){changeInRegion(x, myrange, refgenome, annot)}, mc.cores = threads)
    pred.full <- rbindlist(predict.list)
    pred.full <- prettyTable(pred.full)
  } else {
    predict.list <- mclapply(files.vcf.gz, function(x){changeAll(x, refgenome, annot)}, mc.cores = threads)
    if (length(predict.list)==1){
      pred.full <- predict.list[[1]]
    } else { 
      pred.full <- rbindlist(predict.list)
    }
    
    pred.full <- prettyTable(pred.full)
  }
  
  # make all variants full table
  vcf.table <- vcf2table(files.vcf)
  colnames(vcf.table) <- c("#CHROM", "start", "REF", "ALT", "sample")
  
  # join predict table and vcf (vars) table by sample and position
  #
  a.list <- list()
  
  var.split <- split.data.frame(vcf.table, vcf.table$sample)
  pred.split <- split.data.frame(pred.full, pred.full$sample)
  # WARNING! if there are positions occured multiple times (with different nuc substitutions)
  # left_join by "start" will not work properly
  # MUST be like so - left_join(by=c("start", "ALT"))
  if (length(var.split) == length(pred.split)){
    # joining
    for(i in c(1:length(var.split))){
      a <- left_join(var.split[[i]],  pred.split[[i]], by="start")
      a.list[[i]] <- a
    }
  } else {
    print("Error! Annotation table and SNP table have different number of samples!")
    break()
  }
  
  
  a.full <- rbindlist(a.list)
  a.full$GENEID <- if_else(is.na(a.full$GENEID), "intergenic", a.full$GENEID)
  # a.full <- select(a.full, -c(width, REF.y, ALT.y, sample.y))
  a.full <- a.full[,-c(6, 8, 9, 19)]
  colnames(a.full) <- c("Chrom", "Pos", "Ref", "Alt", "sample", "strand", "CDSLOC.start", "CDSLOC.end", "PROTEINLOC", "GENEID", "CONSEQUENCE", "REFCODON", "VARCODON", "REFAA", "VARAA")
  
  # spread to wide and arrange by position
  a.full.wide <- spread(a.full, key="sample", value="Alt") %>%  arrange(Pos) #%>%  distinct(REF.x, start, GENEID.y, .keep_all=T)%>% arrange(start)
  
  # annotate rRNA genes (the following coords are suitable for M.tuberculosis only)
  if (species=="Mycobacterium tuberculosis"){
    # rrs
    a.full.wide$GENEID <- if_else(a.full.wide$Pos >= 1471846 & a.full.wide$Pos <= 1473382, "rrs", a.full.wide$GENEID)
    # rrL
    a.full.wide$GENEID <- if_else(a.full.wide$Pos >= 1473658 & a.full.wide$Pos <= 1476795, "rrl", a.full.wide$GENEID)
    # rrf
    a.full.wide$GENEID <- if_else(a.full.wide$Pos >= 1476899 & a.full.wide$Pos <= 1477013, "rrf", a.full.wide$GENEID)
  } else {
    print("rRNA genes will not be annotated!")
  }
  
  
  # deal with the cases when there are different amino acids in different samples because of different nuc substitutions in the same position
  
  left.part <- a.full.wide[,1:8] # Chrom, Pos, Ref, Strand
  ncols0 <- length(a.full.wide) - 1 # 0-based counting in data table
  
  # make right part which will handle ambigous AA and NUC changes
  print("dealing with different substitutions...")
  right.part <- setDT(a.full.wide)[, lapply(.SD[,8:ncols0], function(x) paste(x[!is.na(x)], collapse='/')) , Pos]
  
  # bind both parts & remove redundancy
  a.full.wide.final <- full_join(left.part, right.part, by="Pos") %>% distinct(Pos, Ref, CDSLOC.start, .keep_all = TRUE)
  
  write.table(a.full.wide.final, paste(outTable, ".tsv", sep=""), sep="\t", quote=F, na="", row.names=F)
  # if (writeXLSX){
  #   print("Warning! Writing xlsx may take very long time!..")
  #   write.xlsx(a.full.wide.final,  paste(outTable, ".xlsx", sep=""), showNA = F, row.names = F)
  # }
  # pred.full <- spread(pred.full, key="sample", value="ALT") %>% arrange(start)
  # write.table(pred.full, "aa_changes_wide.tsv", sep="\t", row.names = F, quote=F, na="")
  # write.xlsx(pred.full,  "aa_changes_wide.pred.fullspred.full", showNA = F)
  return(list(a.full, a.full.wide.final))
}
