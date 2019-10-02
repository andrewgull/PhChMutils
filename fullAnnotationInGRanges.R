fullAnnotationInGRanges <- function(resistance_table){
  source("/home/agulyaev/bin/makeAnnotation.R")
  source("/home/agulyaev/bin/resist2GRanges.R")
  resist_ranges <- resist2GRanges(resistance_table)
  annots <- makeAnnotation(withRanges=TRUE, inRanges=resist_ranges, outTable="annotation_rpoABC")
  
}
