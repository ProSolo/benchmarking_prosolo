#!/usr/bin/env Rscript

library('tidyverse')
library('argparser')

argv <- arg_parser("Calculate expected value of amplification error from unfiltered probabilities.") %>%
        add_argument("--calls", help="file with calls for one single cell") %>%
        add_argument("--out", help="output file") %>%
        parse_args();

chrom_levels <- c( "chr1", "chr10", "chr11", "chr12", "chr13", "chr14",
                   "chr15", "chr16", "chr17", "chr18", "chr19", "chr2", "chr20",
                   "chr21", "chr22", "chr3", "chr4", "chr5", "chr6", "chr7",
                   "chr8", "chr9", "chrX")
chrom_exclude_levels <- c( "chrY" )


phred2prob <- function(q) {
  return( 10^( (-q)/10 ) )
}

results <- as_tibble(
  read_tsv(
    argv$calls,
    col_types = "ci__ddddddd__",
    col_names = TRUE,
    na = c("." )
    ) ) %>%
  mutate(CHROM = factor(CHROM, levels = chrom_levels, exclude = chrom_exclude_levels), 
         PROB_ERR_REF = phred2prob(PROB_ERR_REF),
         PROB_ERR_ALT = phred2prob(PROB_ERR_ALT),
         ) %>%
  transmute(
         prob_ERR_set = PROB_ERR_REF + PROB_ERR_ALT
         ) %>%
  summarise(err_rate = mean(prob_ERR_set) )

write_tsv(results, argv$out)

