library(tidyverse)
library('RColorBrewer')
library('extrafont')
loadfonts(device = "pdf")


current <- tibble()

file_l_gt <- str_c("Laehnemann2017/prosolo/single-cells.PNG.alt_sites_only.prosolo.min_sc_cov_1.ado-rates-against-ground-truth.tsv")
ado_rates <- as_tibble( read_tsv(file_l_gt) ) %>%
            add_column( dataset = "granulocytes",
                        type = "ProSolo hom at ground truth het"
            )

file_d_gt <- str_c("Dong2017/prosolo/single-cells.Clones.alt_sites_only.prosolo.min_sc_cov_1.minGQ_30.ado-rates-against-ground-truth.tsv")
current <- as_tibble( read_tsv(file_d_gt) ) %>%
            add_column( dataset = "Dong 2017",
                        type = "ProSolo hom at ground truth het"
            )
ado_rates <- bind_rows(ado_rates, current)

file_known <- str_c("known_ados.tsv")
current <- as_tibble( read_tsv(file_known) ) %>%
            mutate( dataset = citation,
                    ado_rate = ADO) %>%
            select( dataset, ado_rate ) %>%
            add_column( cell = "published",
                        type = "published rate"
            )
ado_rates <- bind_rows(ado_rates, current)

cells_l <- c("PAG1", "PAG2", "PAG5", "PAG9", "PAG10")
for (cell in cells_l) {
  file <- str_c("Laehnemann2017/prosolo/ado-rate/", cell, ".PNG.alt_sites_only.prosolo.min_sc_cov_1.ground_truth_hets_only.expected-ado-rate.tsv")
  current <- as_tibble( read_tsv(file) ) %>%
              add_column( dataset = "granulocytes",
                          type = "ProSolo expected value ADO events",
                          cell = cell
              )
  ado_rates <- bind_rows(ado_rates, current)
}
for (cell in cells_l) {
  file <- str_c("Laehnemann2017/", cell, ".CCS_1_W.ground_truth.het_snps_only.allele_dropout_rate.tsv")
  current <- as_tibble( read_tsv(file) ) %>%
              transmute(ado_rate = ALLELE_DROPOUT_RATE) %>%
              add_column( dataset = "granulocytes",
                          type = "no alt/ref read at ground truth het",
                          cell = cell
              )
  ado_rates <- bind_rows(ado_rates, current)
}


cells_d <- c("IL-11", "IL-12")
for (cell in cells_d) {
  file <- str_c("Dong2017/prosolo/ado-rate/", cell, ".Clones.alt_sites_only.prosolo.min_sc_cov_1.minGQ_30.ground_truth_hets_only.expected-ado-rate.tsv")
  current <- as_tibble( read_tsv(file) ) %>%
              add_column( dataset = "Dong 2017",
                          type = "ProSolo expected value ADO events",
                          cell = cell
              )
  ado_rates <- bind_rows(ado_rates, current)
}
for (cell in cells_d) {
  file <- str_c("Dong2017/", cell, ".IL-1c.bulk.ground_truth.minGQ_30.het_snps_only.allele_dropout_rate.tsv")
  current <- as_tibble( read_tsv(file) ) %>%
              transmute(ado_rate = ALLELE_DROPOUT_RATE) %>%
              add_column( dataset = "Dong 2017",
                          type = "no alt/ref read at ground truth het",
                          cell = cell
              )
  ado_rates <- bind_rows(ado_rates, current)
}

ado_rates <- ado_rates %>%
    mutate(cell = factor(cell, ordered = TRUE),
           dataset = factor(dataset, ordered = TRUE,
                        levels = c(
                          "Dong 2017",
                          "granulocytes",
                          "Wang 2014",
                          "Hou 2012",
                          "Xu 2012",
                          "Ling 2009",
                          "Spits 2006",
                          "Renwick 2006",
                          "Lodato 2015"
                          )
                      ),
           type = factor(type, ordered = TRUE,
                      levels = c(
                          "ProSolo hom at ground truth het",
                          "ProSolo expected value ADO events",
                          "no alt/ref read at ground truth het",
                          "published rate"
                          )
                        )
           )

type_pal <- c(
  "ProSolo hom at ground truth het" = gray(0),
  "ProSolo expected value ADO events" = gray(0.4),
  "no alt/ref read at ground truth het" = gray(0.7),
  "published rate" = gray(0.9)
)

dataset_pal <- c(
  "granulocytes" = 11,
  "Dong 2017" = 4,
  "Wang 2014" = 0,
  "Hou 2012" = 1,
  "Xu 2012" = 2,
  "Ling 2009" = 10,
  "Spits 2006" = 8,
  "Renwick 2006" = 5,
  "Lodato 2015" = 6
)

ggplot(ado_rates %>% group_by(dataset), aes(x = cell, y = ado_rate) ) +
  geom_col( aes(fill = type), position = position_dodge2(width = .9, preserve = "total") ) +
  geom_point( aes(shape = dataset), size = 3, position = position_dodge2(width = .9, preserve = "total") )+
  scale_fill_manual(
    values = type_pal
    ) +
  scale_shape_manual(
    values = dataset_pal
    ) +
  scale_y_continuous(
    breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25),
    labels = scales::percent_format(accuracy = 1, suffix = "")
    )  +
  theme_bw(base_size=24, base_family="Lato")  +
  theme( axis.text.x = element_text( angle=45, hjust = 1) ) +
  guides(
         fill = guide_legend("allele dropout calculation")
         ) +
  labs(
       x = "single cell / dataset",
       y = "% allele dropout"
       )
ggsave("prosolo_ado-rate_expected-vs-ground-truth.pdf", device = cairo_pdf, width=16, height=7)

