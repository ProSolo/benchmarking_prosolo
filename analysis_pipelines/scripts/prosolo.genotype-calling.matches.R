#!/usr/bin/env Rscript

library('tidyverse')
library('argparser')

argv <- arg_parser("Genotype matching calculations on ProSolo results FDR filtered by ALT allele presence.") %>%
        add_argument("--gt",
		     help="file with ground truth genotypes (incl. path), contents: `CHROM\tPOS\tGT`") %>%
        add_argument("--calls", help="file with calls for one single cell") %>%
        add_argument("--out", help="output file") %>%
        parse_args();

chrom_levels <- c( "chr1", "chr10", "chr11", "chr12", "chr13", "chr14",
                   "chr15", "chr16", "chr17", "chr18", "chr19", "chr2", "chr20",
                   "chr21", "chr22", "chr3", "chr4", "chr5", "chr6", "chr7",
                   "chr8", "chr9", "chrX")
chrom_exclude_levels <- c( "chrY" )


PED_GTs <- as_tibble( read_tsv( argv$gt, col_names = c("CHROM", "POS", "PED_GT") ) ) %>%
  mutate( PED_GT = factor(PED_GT, levels = c("1/1", "0/1", "0/0", NA) ),
          CHROM = factor(CHROM, levels = chrom_levels, exclude = chrom_exclude_levels)
          )

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
         PROB_HOM_REF = phred2prob(PROB_HOM_REF),
         PROB_ADO_TO_REF = phred2prob(PROB_ADO_TO_REF),
         PROB_ERR_REF = phred2prob(PROB_ERR_REF),
         PROB_HOM_ALT = phred2prob(PROB_HOM_ALT),
         PROB_ADO_TO_ALT = phred2prob(PROB_ADO_TO_ALT),
         PROB_ERR_ALT = phred2prob(PROB_ERR_ALT),
         PROB_HET = phred2prob(PROB_HET),
         prob_HOM_REF_set = PROB_HOM_REF + PROB_ERR_ALT,
         prob_HET_set = PROB_HET + PROB_ADO_TO_ALT + PROB_ADO_TO_REF,
         prob_HOM_ALT_set = PROB_HOM_ALT + PROB_ERR_REF,
         max = pmax( prob_HOM_REF_set, prob_HET_set, prob_HOM_ALT_set ),
         GT = case_when(
            max == prob_HOM_REF_set ~ "0/0",
            max == prob_HET_set ~ "0/1",
            max == prob_HOM_ALT_set ~ "1/1"
            ),
	 GT = factor(GT, levels = c("1/1", "0/1", "0/0", NA), exclude = NULL )
         ) %>%
  select( CHROM, POS, GT ) %>%
  add_count( CHROM, POS ) %>%
  distinct( CHROM, POS, .keep_all = TRUE ) %>%
  mutate( GT = if_else( n == 1, factor(GT) , factor(NA_character_) ) ) %>%
  select( CHROM, POS, GT ) %>%
  right_join(PED_GTs, by = c("CHROM", "POS")) %>%
  count( PED_GT, GT )

write_tsv(results, argv$out)

