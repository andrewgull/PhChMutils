#reformat *.snps files to vcf for annotation
# by default use transformed (`sed  "s/:/\t/g" $f | sed  "s/%//g" > $f"t"`) snps files
snps2vcf <- function(ptrn=".snpst"){
  library(data.table)
  library(dplyr)
  # only.snps.list <- lapply(paths, function(path){
  # if (grepl('snps', path)){
  # g <- regexpr("[S, E]RR[0-9]*.", path)
  # sample.name <- regmatches(path, g)
  # x <- read.table(path, head=T,sep="\t")
  # x$sample <- rep(sample.name, nrow(x))
  # return(x)
  # }
  # })
  
  # the variant above is for multiple different paths to snps files
  # find proper transformed files and filter them
  snps.files <- dir(pattern = ptrn)
  
  filter.snps <- function(f){
    # reads files, filters and writes filtered versions of them
    snps <- read.delim(f)
    snps <- snps[snps$Freq >= 90 & snps$Cov >= 10,]
    write.table(snps, paste0(f, "f", sep=""), quote = F, row.names = F, sep = "\t")
  }
  
  filter.snps <- function(f){
    # reads files, filters and writes filtered versions of them
    snps <- read.delim(f)
    snps <- filter_(na.omit(snps), "Freq">=90 & "Cov" >= 10)
    #snps <- snps[snps$Freq >= 90 & snps$Cov >= 10,]
    #snps <- snps[,c(1,2,3,4)]
    #snps$Chrom <- sub("gi|448814763|ref|", "", snps$Chrom, fixed=T)
    #snps$Chrom <- sub("|", "", snps$Chrom, fixed=T)
    write.table(snps, paste0(f, "f", sep=""), quote = F, row.names = F, sep = "\t")
    #return(snps)
  }
  
  
  sapply(snps.files, filter.snps)
  
  # find proper files (they have to be both transformed and filtered) and make sample names
  new_ptrn <- paste(ptrn, "f", sep="")
  files <- dir(pattern = new_ptrn)
  only.snps.list <- lapply(files, function(f){
    x <- read.table(f, head=T, sep="\t")
    sample.name <- sub(new_ptrn, "", f, fixed=TRUE)
    x$sample <- rep(sample.name, nrow(x))
    return(x)
  })
  
  only.snps.df <- rbindlist(only.snps.list)
  # substitute some text variables
  only.snps.df$sample <- sub(".", "", only.snps.df$sample, fixed=TRUE)
  only.snps.df$Chrom <- sub("gi|448814763|ref|", "", only.snps.df$Chrom, fixed = T)
  only.snps.df$Chrom <- sub(".3|", "", only.snps.df$Chrom, fixed = T)
  # select needed columns and rename them
  only.snps.df <- dplyr::select(only.snps.df, c(Chrom, Position, Ref, Var, SamplesRef, sample))
  colnames(only.snps.df) <- c("#CHROM", "POS", "REF", "ALT", "QUAL", "sample")
  
  ################################
  
  # paste ID column in right place
  only.snps.df$ID <- rep(".", nrow(only.snps.df))
  only.snps.df <- only.snps.df  %>% dplyr::select(ID, everything())
  only.snps.df <- only.snps.df  %>% dplyr::select(c(`#CHROM`, POS), everything())
  # add more column
  only.snps.df$FILTER <- rep("0", nrow(only.snps.df))
  only.snps.df$INFO <- rep("AB=0;ABP=0;", nrow(only.snps.df))
  only.snps.df$FORMAT <- rep("GT:DP:AD", nrow(only.snps.df))
  only.snps.df$unknown <- rep("1:126:0", nrow(only.snps.df))
  
  
  # write to files
  only.snps.list <- split(only.snps.df, only.snps.df$sample)
  #x <- lapply(df.list, function(item){write.table(item, paste(item[1,7], "vcf", sep="."), quote=F, row.names = F, sep="\t")})
  x <- lapply(only.snps.list, function(item){
    id <- item$sample[1]
    df <- dplyr::select(item, -sample)
    filename <- paste(id, "vcf", sep=".")
    write.table(df, filename, quote=F, row.names = F, sep="\t")
  })
  #rm(x)
}
