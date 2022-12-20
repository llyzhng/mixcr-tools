# Sample UMIs from umi_stats.txt and output in formatted tag_table.tsv for MiTool
library(tidyverse)
library(optparse)

# ===================================================================================== #
# Command line arguments

option_list = list(
  make_option(c("-i", "--umi_stats"), type="character", default=NULL,
              help="umi-stats.txt from MiTool", metavar="character"),
  make_option(c("-o", "--tag_table"), type="character", default=NULL,
              help="tag-table.tsv for MiTool", metavar="character"),
  make_option(c("-p", "--sampled_umi_stats"), type="character", default=NULL,
              help="umi-stats.txt for sampled UMIs. Not output if not supplied.", metavar="character"),
  make_option(c("-f", "--filtered_mic"), type="character", default="filtered.mic",
              help="path for filtered.mic containing all the reads with selected UMIs", metavar="character"),
  make_option(c("-s", "--seed"), type="numeric", default=NULL,
              help="random seed", metavar="numeric"),
  make_option(c("-n", "--num_to_sample"), type="numeric", default=NULL,
              help="number of UMIs to sample", metavar="numeric")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# ===================================================================================== #

umis <- read.delim(opt$umi_stats) %>%
  as_tibble()

set.seed(opt$seed)
sampled_nums <- sample(1:nrow(umis), size = opt$num_to_sample, replace = F)
sampled_umis <- umis[sampled_nums, ]

if (!is.null(opt$sampled_umi_stats)) {
  write.table(sampled_umis, opt$sampled_umi_stats, sep = "\t", quote = F, col.names = T, row.names = F)
}

tag_table <- tibble(UMI = c("UMI", sampled_umis$UMI), 
                    file = c("", rep(opt$filtered_mic, opt$num_to_sample)))
write.table(tag_table, opt$tag_table, sep = "\t", quote = F, col.names = F, row.names = F)