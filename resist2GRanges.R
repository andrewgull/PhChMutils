resist2GRanges <- function(resdf_name, row_max=1370){
  # a function to make IRanges object of AB-resistance genes 
  # first three column names have to be: "Genome position", "Gene", "AA exchange"
  
  # inRanges - GRanges object of regions of interest
  # example: GRanges(seqnames = "NC_000962", ranges = IRanges(start=c(3877464, 759807, 763370), end=c(3878507, 763325, 767320), names=c("rpoA", "rpoB", "rpoC")))

  library(GenomicRanges)
  library(IRanges)
  #library(readxl)
  
  resdf <- read.csv(resdf_name, nrows = row_max, sep='\t')
  resdf <- resdf[,c(1:3)]
  genes <- unique(resdf[,2])
  res.list <- split(resdf, resdf$Gene)
  gene.starts <- sapply(res.list, function(x){min(x[,1])})
  gene.ends <- sapply(res.list, function(x){max(x[,1])})
  range_obj <- GRanges(seqnames="NC_000962", ranges=IRanges(start=gene.starts, end=gene.ends, names=sort(genes)))
  return(range_obj)
}