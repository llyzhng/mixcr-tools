library(optparse)
library(tidyverse)

# ===================================================================================== #
# Command line arguments

option_list = list(
  make_option(c("-x", "--table1"), type="character", default=NULL,
              help="first clone table to combine", metavar="character"),
  make_option(c("-y", "--table2"), type="character", default=NULL,
              help="second clone table to combine", metavar="character"),
  make_option(c("-o", "--outfile"), type="numeric", default=NULL,
              help="path for output combined tabled", metavar="numeric")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# ===================================================================================== #

t1 <- read.delim(opt$table1)
t2 <- read.delim(opt$table2)

comb <- rbind(t1, t2)
totalReads <- sum(comb$readCount)
comb$readFraction <- comb$readCount / totalReads

if (any(grepl("uniqueUMI", colnames(comb)))) {
  totalUMI <- sum(comb$uniqueUMICount)
  comb$uniqueUMIFraction <- comb$uniqueUMICount / totalUMI
} 

write.table(comb, opt$outfile, 
            sep = "\t", quote = F, col.names = T, row.names = F)
