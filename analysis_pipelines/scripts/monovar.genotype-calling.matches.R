#!/usr/bin/env Rscript

library('tidyverse')
library('argparser')

argv <- arg_parser("Genotype matching calculations on MonoVar.") %>%
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

results <- as_tibble(
  read_tsv(
    argv$calls,
    skip = 1,
    col_types = "ci__c",
    col_names = c( "CHROM", "POS", "GT" ) ) ) %>%
  mutate(
    CHROM = factor(CHROM, levels = chrom_levels, exclude = chrom_exclude_levels),
    GT = factor(GT, levels = c("1/1", "0/1", "0/0", NA) ) ) %>%
  right_join(PED_GTs, by = c("CHROM", "POS")) %>%
  count( PED_GT, GT )

write_tsv(results, argv$out)

