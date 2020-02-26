#!/usr/bin/env Rscript

library('tidyverse')
library('argparser')

argv <- arg_parser("Plot empirical cumulative density of amplification error probabilities of prosolo calls.") %>%
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

err_probs <- as_tibble(
  read_tsv(
    argv$calls,
    col_types = "ci__ddddddd__",
    col_names = TRUE,
    na = c("." )
    ) ) %>%
  mutate(CHROM = factor(CHROM, levels = chrom_levels, exclude = chrom_exclude_levels), 
         PROB_HOM_REF = phred2prob(PROB_HOM_REF),
         PROB_ADO_TO_REF = phred2prob(PROB_ADO_TO_REF),
         PROB_ERR_REF = phred2prob(PROB_ERR_REF),
         PROB_HOM_ALT = phred2prob(PROB_HOM_ALT),
         PROB_ADO_TO_ALT = phred2prob(PROB_ADO_TO_ALT),
         PROB_ERR_ALT = phred2prob(PROB_ERR_ALT),
         PROB_HET = phred2prob(PROB_HET)
         ) %>%
  transmute(
         prob_ERR_set = PROB_ERR_REF + PROB_ERR_ALT
         )

ggplot(err_probs, aes(prob_ERR_set ) ) +
  stat_ecdf() +
  theme_bw()
ggsave(argv$out, device = svg, width = 9, height = 6)

