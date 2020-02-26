#!/usr/bin/env Rscript

library('tidyverse')
library('argparser')

argv <- arg_parser("True and false positive / negative calculations on ProSolo results FDR filtered by ALT allele presence.") %>%
        add_argument("--gt",
		     help="file with ground truth genotypes (incl. path), contents: `CHROM\tPOS\tGT`") %>%
        add_argument("--calls", help="file with calls for one single cell") %>%
        add_argument("--out", help="output file") %>%
        parse_args();

PED_GTs <- as_tibble( read_tsv( argv$gt, col_names = c("CHROM", "POS", "PED_GT") ) ) %>%
  mutate( PED_GT = factor(PED_GT, levels = c("1/1", "0/1", "0/0", NA) ),
          CHROM = factor(CHROM) )

phred2prob <- function(q) {
  return( 10^( (-q)/10 ) )
}

probs <- as_tibble(
  read_tsv(
    argv$calls,
    col_types = "ci__ddddddd__",
    col_names = TRUE,
    na = c("." )
    ) ) %>%
  mutate(CHROM = factor(CHROM), 
         PROB_HOM_REF = phred2prob(PROB_HOM_REF),
         PROB_ADO_TO_REF = phred2prob(PROB_ADO_TO_REF),
         PROB_ERR_REF = phred2prob(PROB_ERR_REF),
         PROB_HOM_ALT = phred2prob(PROB_HOM_ALT),
         PROB_ADO_TO_ALT = phred2prob(PROB_ADO_TO_ALT),
         PROB_ERR_ALT = phred2prob(PROB_ERR_ALT),
         PROB_HET = phred2prob(PROB_HET),
         prob_ref_set = PROB_HOM_REF + PROB_ERR_ALT,
         prob_alt_set = PROB_HOM_ALT + PROB_ERR_REF + PROB_HET + PROB_ADO_TO_ALT + PROB_ADO_TO_REF,
         call = if_else( prob_ref_set < prob_alt_set, "ALT", "REF"),
	 call = factor(call, levels = c("REF", "ALT") ) ) %>%
  select( CHROM, POS, call ) %>%
  group_by( CHROM, POS ) %>%
  summarise( call = if_else( any( call == "ALT" ), "ALT", "REF" ) ) %>%
  ungroup() %>%
  right_join(PED_GTs, by = c("CHROM", "POS")) %>%
  mutate(ped_alt_set_match = case_when(
           is.na( PED_GT ) ~ "no_ground_truth",
           is.na( call ) & (PED_GT == "0/1" | PED_GT == "1/1") ~ "P_not_called",
           is.na( call ) & (PED_GT == "0/0") ~ "N_not_called",
           call == "ALT" & (PED_GT == "0/1" | PED_GT == "1/1") ~ "TP",
           call == "ALT" & (PED_GT == "0/0") ~ "FP",
           call == "REF" & (PED_GT == "0/1" | PED_GT == "1/1") ~ "FN",
           call == "REF" & (PED_GT == "0/0") ~ "TN" ),
         ped_alt_set_match = 
           factor( ped_alt_set_match,
                   levels=c("TP", "FP", "FN", "TN", "P_not_called","N_not_called", "no_ground_truth") )
	 ) %>%
  select( CHROM, POS, ped_alt_set_match )

## ALT presence calculations

mls_pn <- probs %>%
  count(ped_alt_set_match) %>%
  complete(ped_alt_set_match) %>%
  mutate(n=ifelse(is.na(n),0,n)) %>%
  spread(ped_alt_set_match, n )

write_tsv(mls_pn, argv$out)

