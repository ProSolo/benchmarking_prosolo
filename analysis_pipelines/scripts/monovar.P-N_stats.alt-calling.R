#!/usr/bin/env Rscript

library('tidyverse')
library('argparser')

argv <- arg_parser("True and false positive / negative calculations on MonoVar.") %>%
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
    col_types = "ci__c",
    col_names = c( "CHROM", "POS", "GT" ) ) ) %>%
  mutate(
    CHROM = factor(CHROM), 
    GT = factor(GT, levels = c("1/1", "0/1", "0/0", NA) ) ) %>%
  right_join(PED_GTs, by = c("CHROM", "POS")) %>%
  mutate(
    monovar_alt_match = case_when(
      is.na( PED_GT ) ~ "no_ground_truth",
      is.na( GT ) & (PED_GT == "0/1" | PED_GT == "1/1") ~ "P_not_called",
      is.na( GT ) & (PED_GT == "0/0") ~ "N_not_called",
      (GT == "0/1" | GT == "1/1") & (PED_GT == "0/1" | PED_GT == "1/1") ~ "TP",
      (GT == "0/1" | GT == "1/1") & (PED_GT == "0/0") ~ "FP",
      GT == "0/0" & (PED_GT == "0/1" | PED_GT == "1/1") ~ "FN",
      GT == "0/0" & (PED_GT == "0/0") ~ "TN"),
    monovar_alt_match = factor( monovar_alt_match,
                   levels=c("TP", "FP", "FN", "TN", "P_not_called","N_not_called", "no_ground_truth") ) )

## ALT presence calculations

pn <- annotated %>%
  count(monovar_alt_match) %>%
  complete(monovar_alt_match) %>%
  mutate(n=ifelse(is.na(n),0,n)) %>%
  spread(monovar_alt_match, n )

write_tsv(pn, argv$out)

