library(tidyverse)
library(optparse)

# ===================================================================================== #
# Command line arguments

option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL,
              help="path for clone table to collapse", metavar="character"),
  make_option(c("-o", "--out_heavy"), type="character", default=NULL,
              help="path for csv with output heavy chain clones", metavar="character"),
  make_option(c("-l", "--out_light"), type="character", default=NULL,
              help="path for csv with output light chain clones to keep clone", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# ===================================================================================== #
# Functions

getName <- function(x) {
  # remove allele
  gsub("\\*.*", "", x)
}

generateConsensusClones <- function(clone_tb) {
  # Remove allele information from V and J segment calls 
  clone_sub <- clone_tb %>%
    select(v_call, j_call, cdr3, clone_id) %>%
    mutate(bestVGene = getName(v_call),
           bestJGene = getName(j_call))
  
  # Group all CDR3 information (sequence and UMI counts) for clones
  clone_cdr3 <- clone_sub %>%
    select(clone_id, bestVGene, bestJGene, cdr3) %>%
    group_by(clone_id, bestVGene, bestJGene, cdr3) %>%
    mutate(cdr3_umi_count = n()) %>%
    distinct() %>%
    arrange(desc(cdr3_umi_count)) 
  best_cdr3 <- clone_cdr3 %>%
    group_by(clone_id) %>%
    filter(cdr3_umi_count == max(cdr3_umi_count)) %>%
    rename(bestCDR3 = cdr3) %>%
    ungroup() %>%
    select(clone_id, bestVGene, bestJGene, bestCDR3)
  all_cdr3 <- clone_cdr3 %>%
    group_by(clone_id, bestVGene, bestJGene) %>%
    mutate(cdr3_info = paste0(cdr3, ":", cdr3_umi_count)) %>%
    mutate(all_cdr3 = paste(cdr3_info, collapse = ",")) %>%
    group_by(clone_id, bestVGene, bestJGene, all_cdr3) %>%
    summarize(uniqueUMICount = sum(cdr3_umi_count)) %>%
    ungroup() %>%
    arrange(desc(uniqueUMICount))
  comb_cdr3 <- left_join(best_cdr3, all_cdr3)
  
  # Summarize clones with best V/J/CDR3, all CDR3 info, UMI count and fraction
  consensus_clones <- comb_cdr3 %>%
    mutate(uniqueUMIFraction = uniqueUMICount / sum(uniqueUMICount)) %>%
    arrange(desc(uniqueUMIFraction)) %>% 
    select(bestVGene, bestJGene, bestCDR3, all_cdr3, uniqueUMICount, uniqueUMIFraction)
  
  return(consensus_clones)
}
# ===================================================================================== #

clone_pass <- read.delim(opt$infile) %>%
  as_tibble() 

clones_heavy <- clone_pass %>%
  filter(locus == "IGH") 

if (nrow(clones_heavy) >= 1) {
  clones_heavy_consensus <- clones_heavy %>% 
  generateConsensusClones()
  write.table(clones_heavy_consensus, opt$out_heavy,
              sep = "\t", quote = F, row.names = F, col.names = T)
} else {
  message("No heavy chain clones")
}

clones_light <- clone_pass %>%
  filter(locus == "IGL" | locus == "IGK") 

if (nrow(clones_light) >= 1) {
  clones_light_consensus <- clones_light %>% 
    generateConsensusClones()
  write.table(clones_light_consensus, opt$out_light,
              sep = "\t", quote = F, row.names = F, col.names = T)
} else {
  message("No light chain clones")
}

