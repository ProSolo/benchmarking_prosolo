library('tidyverse')
library('stringr')
library('ggthemes')
library('scales')
library('RColorBrewer')

Dong_SCs <- c("IL-11","IL-12")

remove(current)
remove(counts)

# read in ProSolo data

fdrs <- c( "0-2" )

for ( sc in Dong_SCs ) {
  for (f in fdrs ) {
    file <- str_c("prosolo/fdr_alt-presence/", sc,".Clones.alt_sites_only.fdr_", f , "_alt-presence.prosolo.min_sc_cov_1.minGQ_30.gt_matching.ML_set.tsv")
    current <- as_tibble( read_tsv( file ) ) %>%
                add_column( cell = sc, 
                            software = "ProSolo, default",
                            filter = if_else(f == "1", 1.0, as.numeric( str_replace(f, '-', '.') )) )
    if (!exists("counts")) {
      counts <- current
    } else {
      counts <- bind_rows(counts,current)
    }
  }
}

# Read in Monovar data
for ( sc in Dong_SCs ) {
  file <- str_c("monovar/", sc,".single_cells.1.c0.monovar.minGQ_30.gt_matching.tsv")
  current <- as_tibble( read_tsv( file ) ) %>%
              add_column(cell = sc,
                         software = "MonoVar, no cons",
                         filter = 1)
  if (!exists("counts")) {
    counts <- current
  } else {
    counts <- bind_rows(counts,current)
  }
}

GT_levels <- c("not called", "1/1", "0/1", "0/0")

# order factor levels of GT and PED_GT for consistent bar stacking
counts <- counts %>%
  mutate( match = case_when(
                          is.na(GT) ~ "not called",
                          GT == PED_GT ~ "correct GT",
                          TRUE ~ "false GT"
                 ),
          match = factor(match, levels = c("not called", "false GT", "correct GT" ) ),
          GT = replace_na(GT, "not called"),
          PED_GT = replace_na(PED_GT, "not called"),
          PED_GT = factor( PED_GT, levels = GT_levels, exclude = NULL, ordered = TRUE),
          GT = factor( GT, levels = GT_levels, exclude = NULL, ordered = TRUE )
  )

# color naming for colorblind palette
Black <- colorblind_pal()(8)[[1]]
Orange <- colorblind_pal()(8)[[2]]
Sky_Blue <- colorblind_pal()(8)[[3]]
bluish_Green <- colorblind_pal()(8)[[4]]
Yellow <- colorblind_pal()(8)[[5]]
Blue <- colorblind_pal()(8)[[6]]
Vermillion_Red <- colorblind_pal()(8)[[7]]
reddish_Purple <- colorblind_pal()(8)[[8]]

# color assignment to different sets
hom_ref <- reddish_Purple
het <- Sky_Blue
hom_alt <- Blue
error <- Vermillion_Red
correct <- bluish_Green
not_called <- Yellow

true_false_no_call_palette <- c( not_called, error, correct)
true_false_palette <- c( error, correct)

gt_palette <- c(
  "1/1" = hom_alt,
  "0/1" = het,
  "0/0" = hom_ref,
  "not called" = not_called )

ggplot(counts %>% filter(PED_GT != "0/0" & PED_GT != "not called"), aes( x = interaction(software, filter, sep = ", "), y = n )) +
  facet_grid(PED_GT ~ cell) +
  geom_col( aes(
              #colour = GT,
              fill = match
              ),
            width = 0.6,
            size = 1.2 ) + 
  theme_bw( ) +
  theme( axis.text.x = element_text( angle=45, hjust = 1) ) +
  labs( title = "recall of pedigree ground truth genotypes",
        x = "software and filter level",
        y = "number of calls",
        tag = "A"
        ) +
  guides( 
          #colour = guide_legend("caller genotype", override.aes = list(fill = "white"), order = 1),
          fill = guide_legend("ground truth matching", order = 2)
         ) +
  #scale_colour_manual(values = gt_palette ) +
  scale_fill_manual(values = true_false_no_call_palette ) +
  scale_x_discrete(limits = c("ProSolo, default, 0.2", "MonoVar, no cons, 1")) +
  scale_y_continuous(labels = comma)
ggsave("Dong2017_genotyping_recall_prosolo_monovar.pdf", width=5, height=5)

ggplot(counts %>% filter(GT != "not called" & GT != "0/0"), aes( x = interaction(software, filter, sep = ", "), y = n )) +
  facet_grid(GT ~ cell) +
  geom_col( aes(
              #colour = PED_GT,
              fill = match
              ),
            width = 0.6,
            size = 1.2 ) + 
  theme_bw( ) +
  theme( axis.text.x = element_text( angle=45, hjust = 1) ) +
  labs( title = "precision of called genotypes",
        x = "software and filter level",
        y = "number of calls",
        tag = "B"
        ) +
  guides( 
          #colour = guide_legend("ground truth genotype", override.aes = list(fill = "white"), order = 1),
          fill = guide_legend("ground truth matching", order = 2)
         ) +
  #scale_colour_manual(values = gt_palette ) +
  scale_fill_manual(values = true_false_palette ) +
  scale_x_discrete(limits = c("ProSolo, default, 0.2", "MonoVar, no cons, 1")) +
  scale_y_continuous(labels = comma)
ggsave("Dong2017_genotyping_precision_prosolo_monovar.pdf", width=5, height=5)

ggplot(counts %>% filter(PED_GT != "0/0" & PED_GT != "not called"), aes( x = interaction(software, filter, sep = ", "), y = n )) +
  facet_grid(PED_GT ~ cell) +
  geom_col( aes(
              fill = GT
              ),
            width = 0.6,
            size = 1.2 ) + 
  theme_bw( ) +
  theme( axis.text.x = element_text( angle=45, hjust = 1) ) +
  labs( title = "called genotypes for ground truth genotypes",
        x = "software and filter level",
        y = "number of calls",
        tag = "A"
        ) +
  guides( 
          fill = guide_legend("caller genotype", order = 1)
         ) +
  scale_fill_manual(values =  gt_palette ) + #c(gt_palette[[4]], gt_palette[[1]], gt_palette[[2]] ) ) +
  scale_x_discrete(limits = c("ProSolo, default, 0.2", "MonoVar, no cons, 1")) +
  scale_y_continuous(labels = comma)
ggsave("Dong2017_genotyping_calls_per_ground_truth_prosolo_monovar.pdf", width=5, height=5)

ggplot(counts %>% filter(GT != "not called" & GT != "0/0"), aes( x = interaction(software, filter, sep = ", "), y = n )) +
  facet_grid(GT ~ cell) +
  geom_col( aes(
              fill = PED_GT
              ),
            width = 0.6,
            size = 1.2 ) + 
  theme_bw( ) +
  theme( axis.text.x = element_text( angle=45, hjust = 1) ) +
  labs( title = "ground truth for called genotypes",
        x = "software and filter level",
        y = "number of calls",
        tag = "B"
        ) +
  guides( 
          fill = guide_legend("ground truth genotype", order = 1)
         ) +
  scale_fill_manual(values = gt_palette ) +
  scale_x_discrete(limits = c("ProSolo, default, 0.2", "MonoVar, no cons, 1")) +
  scale_y_continuous(labels = comma)
ggsave("Dong2017_genotyping_ground_truth_per_calls_prosolo_monovar.pdf", width=5, height=5)
