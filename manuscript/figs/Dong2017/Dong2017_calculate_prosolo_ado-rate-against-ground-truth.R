library('tidyverse')
library('stringr')

SCs <- c("IL-11","IL-12")

remove(current)
remove(counts)

# read in ProSolo data

for ( sc in SCs ) {
  #for (f in fdrs ) {
    file <- str_c("prosolo/", sc,".Clones.alt_sites_only.prosolo.min_sc_cov_1.minGQ_30.gt_matching.ML_set.tsv")
    current <- as_tibble( read_tsv( file ) ) %>%
                add_column( cell = sc, 
                            software = "ProSolo"#,
                            #filter = if_else(f == "1", 1.0, as.numeric( str_replace(f, '-', '.') ))
                            )
    if (!exists("counts")) {
      counts <- current
    } else {
      counts <- bind_rows(counts,current)
    }
  #}
}

GT_levels <- c("not called", "1/1", "0/1", "0/0")

# order factor levels of GT and PED_GT for consistent bar stacking
total <- counts %>%
  filter(PED_GT == "0/1" & !is.na(GT)) %>%
  group_by(cell) %>%
  summarise(total = sum(n))

ado <- counts %>%
  filter(PED_GT == "0/1" & GT != "0/1" & !is.na(GT)) %>%
  group_by(cell) %>%
  summarise(ado = sum(n))

ado_rates <- full_join(ado, total, by = "cell") %>%
  mutate(ado_rate = ado/total) %>%
  select(cell, ado_rate)
  
write_tsv(ado_rates, "prosolo/single-cells.Clones.alt_sites_only.prosolo.min_sc_cov_1.minGQ_30.ado-rates-against-ground-truth.tsv")

