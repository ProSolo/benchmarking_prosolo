#!/usr/bin/env Rscript

library('tidyverse')
library('argparser')

argv <- arg_parser("True and false positive / negative calculations on SCcaller calls.") %>%
        add_argument("--gt",
		     help="file with ground truth genotypes (incl. path), contents: `CHROM\tPOS\tGT`") %>%
        add_argument("--calls", help="file with calls for one single cell") %>%
        add_argument("--out", help="output file") %>%
        parse_args();

PED_GTs <- as_tibble( read_tsv( argv$gt, col_names = c("CHROM", "POS", "PED_GT") ) ) %>%
  mutate( PED_GT = factor(PED_GT, levels = c("1/1", "0/1", "0/0", NA) ),
          CHROM = factor(CHROM) )

annotated <- as_tibble(
  read_tsv(
    argv$calls,
    skip = 1,
    col_types = "ci___________",
    col_names = c("CHROM", "POS") ) ) %>%
  mutate( CHROM = factor(CHROM) ) %>%
  add_column( call = "ALT" ) %>%
  right_join(PED_GTs, by = c("CHROM", "POS")) %>%
  mutate(
    sccaller_alt_match = case_when(
      is.na( PED_GT ) ~ "no_ground_truth",
      is.na( call ) & (PED_GT == "0/1" | PED_GT == "1/1") ~ "P_not_called",
      is.na( call ) & (PED_GT == "0/0") ~ "N_not_called",
      call == "ALT" & (PED_GT == "0/1" | PED_GT == "1/1") ~ "TP",
      call == "ALT" & (PED_GT == "0/0") ~ "FP"),
    sccaller_alt_match = factor( sccaller_alt_match,
                   levels=c("TP", "FP", "FN", "TN", "P_not_called","N_not_called", "no_ground_truth") ) )

## ALT presence calculations

pn <- annotated %>%
  count(sccaller_alt_match) %>%
  complete(sccaller_alt_match) %>%
  mutate(n=ifelse(is.na(n),0,n)) %>%
  spread(sccaller_alt_match, n )

write_tsv(pn, argv$out)

