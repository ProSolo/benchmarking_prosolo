library('tidyverse')
library('argparser')

argv <- arg_parser("True and false positive / negative calculations on SCAN-SNV data.") %>%
        add_argument("--gt",
		     help="file with ground truth genotypes (incl. path), contents: `CHROM\tPOS\tGT`") %>%
        add_argument("--calls", help="RDA file with calls for one single cell") %>%
        add_argument("--out", help="output file") %>%
        parse_args();

PED_GTs <- as_tibble( read_tsv( argv$gt, col_names = c("CHROM", "POS", "PED_GT") ) ) %>%
  mutate( PED_GT = factor(PED_GT, levels = c("1/1", "0/1", "0/0", NA) ),
          CHROM = factor(CHROM) )

load(argv$calls)

cell <- str_match(argv$calls, "^[^.]+") %>% str_replace('-','.')

annotated <- somatic %>%
  filter( pass == TRUE ) %>%
  rename( CHROM = chr,
          POS = pos) %>%
  mutate(
    CHROM = factor(CHROM), 
    !! cell := factor(!! cell, levels = c("1/1", "0/1", "0/0", NA) ) ) %>%
  right_join(PED_GTs, by = c("CHROM", "POS")) %>%
  mutate(
    scansnv_alt_match = case_when(
      is.na( PED_GT ) ~ "no_ground_truth",
      is.na( !! cell ) & (PED_GT == "0/1" | PED_GT == "1/1") ~ "P_not_called",
      is.na( !! cell ) & (PED_GT == "0/0") ~ "N_not_called",
      (!! cell == "0/1" | !! cell == "1/1") & (PED_GT == "0/1" | PED_GT == "1/1") ~ "TP",
      (!! cell == "0/1" | !! cell == "1/1") & (PED_GT == "0/0") ~ "FP",
      !! cell == "0/0" & (PED_GT == "0/1" | PED_GT == "1/1") ~ "FN",
      !! cell == "0/0" & (PED_GT == "0/0") ~ "TN"),
    scansnv_alt_match = factor( scansnv_alt_match,
                   levels=c("TP", "FP", "FN", "TN", "P_not_called","N_not_called", "no_ground_truth") ) )

## ALT presence calculations

pn <- annotated %>%
  count(scansnv_alt_match) %>%
  complete(scansnv_alt_match) %>%
  mutate(n=ifelse(is.na(n),0,n)) %>%
  spread(scansnv_alt_match, n )

write_tsv(pn, argv$out)
