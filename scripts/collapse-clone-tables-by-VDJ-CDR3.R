library(optparse)

# ===================================================================================== #
# Command line arguments

option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL,
              help="path for clone table to collapse", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="path for output combined tabled", metavar="character"),
  make_option(c("-u", "--umi_threshold"), type="numeric", default=0,
              help="minimum uniqueUMICount to keep clone", metavar="numeric"),
  make_option(c("-r", "--read_threshold"), type="numeric", default=0,
              help="minimum readCount to keep clone", metavar="numeric")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# ===================================================================================== #

x <- read.delim(opt$infile)

# check if table has columns bestVGene, bestJGene, bestDGene
hasGeneCols <- any(any(grepl("bestVGene", colnames(x))),
                 any(grepl("bestJGene", colnames(x))), 
                 any(grepl("bestDGene", colnames(x))))

if (!hasGeneCols) {
  getName <- function(x) {
    gsub("\\*.*", "", x)
  }
  x$bestVGene <- getName(x$allVHitsWithScore)
  x$bestDGene <- getName(x$allDHitsWithScore)
  x$bestJGene <- getName(x$allJHitsWithScore)
}

x$VDJ <- paste(x$bestVGene, x$bestDGene, x$bestJGene, sep = "|")
x$nSeqCDR3[is.na(x$nSeqCDR3)] <- ""

 
if (any(grepl("uniqueUMI", colnames(x)))) {
  # use UMI if available
  collapsed_x <- aggregate(cbind(uniqueUMICount, uniqueUMIFraction) ~ VDJ + nSeqCDR3, data = x, FUN=sum)
  collapsed_x <- collapsed_x[uniqueUMICount >= opt$umi_threshold, ]
  new_sum <- sum(collapsed_x$uniqueUMICount)
  collapsed_x$cloneFraction <- collapsed_x$uniqueUMICount / new_sum
} else {
  # use read count
  collapsed_x <- aggregate(cbind(readCount, readFraction) ~ VDJ + nSeqCDR3, data = x, FUN=sum)
  collapsed_x <- collapsed_x[readCount >= opt$read_threshold, ]
  new_sum <- sum(collapsed_x$readCount)
  collapsed_x$cloneFraction <- collapsed_x$readCount / new_sum
}
colnames(collapsed_x)[4] <- "cloneFraction"

collapsed_x <- collapsed_x[order(-collapsed_x$cloneFraction), ]

write.table(collapsed_x, opt$outfile, 
            sep = "\t", quote = F, row.names = F, col.names = T)
