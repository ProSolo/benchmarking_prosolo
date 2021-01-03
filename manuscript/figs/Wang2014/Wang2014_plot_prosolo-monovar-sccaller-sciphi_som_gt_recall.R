# packages, colors, palettes and stuff
source("../common.R")

tumor_cells <- c(
    "a1",
    "a2",
    "a3",
    "a4",
    "a5",
    "a6",
    "a7",
    "a8",
    "h1",
    "h2",
    "h3",
    "h4",
    "h5",
    "h6",
    "h7",
    "h8"
  )

fdr <- c( "1" )

current <- tibble()
recall_per_cell <- tibble()

# read in ProSolo data
for ( sc in tumor_cells ) {
  for (f in fdr ) {
    file <- str_c(
      "prosolo/fdr_alt-presence/som_gt_recall/",
      sc,
      ".TNBC-Pop-Tumor.alt_sites_only.fdr_",
      f,
      "_alt-presence.prosolo.min_sc_cov_1.gt_calls.ML_set.somatic_genotype_recall.tsv"
    )
    current <- as_tibble(
        read_tsv(
          file,
          col_types = "cicccccc"
          )
      ) %>%
      add_column(
        cell = sc, 
        software = "ProSolo default",
        filter = if_else(f == "1", 1.0, as.numeric( str_replace(f, '-', '.') ))
      )
    recall_per_cell <- bind_rows(recall_per_cell,current)
  }
}

# read in SCIPhI data
for ( sc in tumor_cells ) {
  for (f in fdr ) {
    file <- str_c(
      "sciphi/fdr_alt-presence/som_gt_recall/",
      sc,
      ".tumor_cells.80000-40000_iterations.default.sciphi.alt_prob.fdr_",
      f,
      "_alt-presence.gt_calls.ML_set.somatic_genotype_recall.tsv"
    )
    current <- as_tibble(
        read_tsv(
          file,
          col_types = "cicccccc"
          )
      ) %>%
      add_column(
        cell = sc, 
        software = "SCIPhI default",
        filter = if_else(f == "1", 1.0, as.numeric( str_replace(f, '-', '.') ))
      )
    recall_per_cell <- bind_rows(recall_per_cell,current)
  }
}

# Read in Monovar data
for ( sc in tumor_cells ) {
  for (f in fdr ) {
    file <- str_c(
      "monovar/som_gt_recall/",
      sc,
      ".tumor_cells.",
      f,
      ".c1.monovar.somatic_genotype_recall.tsv"
    )
    current <- as_tibble(
        read_tsv(
          file,
          col_types = "cicccccc"
          )
      ) %>%
      add_column(cell = sc,
        software = "MonoVar default",
        filter = if_else(f == "1", 1.0, as.numeric( str_replace(f, '-', '.') ))
      )
    recall_per_cell <- bind_rows(recall_per_cell,current)
  }
}

# Read in SCcaller data
for ( sc in tumor_cells ) {
  for (f in c("1-0") ) {
    file <- str_c(
      "sccaller/som_gt_recall/",
      sc,
      ".dbsnp.low-cov.sccaller.varcall.cutoff.fil-amp-err_alpha_",
      f,
      ".gt_calls.somatic_genotype_recall.tsv"
    )
    current <- as_tibble(
        read_tsv(
          file,
          col_types = "cicccccc"
          )
      ) %>%
      add_column(
        cell = sc,
        software = "SCcaller sensitive dbsnp",
        filter = if_else(f == "1", 1.0, as.numeric( str_replace(f, '-', '.') ))
      )
    recall_per_cell <- bind_rows(recall_per_cell,current)
  }
}

# Read in validated somatics annotations

validated_somatics <- read_tsv(
  "Wang2014_ground_truth_non_synonymous_variants.hg18_to_hg19.tsv"
)

max_subclonal_het <- validated_somatics %>%
  filter(class == "subclonal") %>%
  filter(zygosity == "het") %>%
  mutate(
    Duplex_Freq = as.double(str_replace(Duplex_Freq, ',', '.'))
  ) %>%
  summarize(max(Duplex_Freq))
write_tsv(max_subclonal_het, path = "Wang2014_max_subclonal_duplex_frequency.tsv")  

expected_numbers <- validated_somatics %>%
  select(class, chrom, pos, REF, VAR, cells, zygosity) %>%
  rename(
    CHROM = chrom,
    POS = pos,
    REF = REF,
    ALT = VAR,
    number_of_cells = cells,
    CLONALITY = class,
    GT = zygosity
  ) %>%
  mutate(
    number_of_cells = if_else(number_of_cells == "NaN", 16, number_of_cells),
    GT = case_when(
      GT == "het" ~ "0/1",
      GT == "hom" ~ "1/1",
      TRUE ~ "unknown"
    )
  ) %>%
  add_column(
    software = "expected"
  )
  
per_software <- recall_per_cell %>%
  filter(GT_CALL == "TP") %>%
  count(software, CHROM, POS, REF, ALT, CLONALITY, GT, name = "number_of_cells") %>%
  full_join(expected_numbers) %>%
  complete(software, nesting(CHROM, POS, REF, ALT, CLONALITY, GT), fill = list(number_of_cells = 0))
  
q <- ggplot(
    per_software
  ) +
  geom_bar(
    aes(
      x = number_of_cells,
      fill = software
    ),
    position = position_dodge( preserve = "single" )
  ) + 
  facet_grid(
    cols = vars(GT),
    rows = vars(CLONALITY),
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
  ) +
  theme_bw( ) +
  labs( title = "somatic genotype recall in TNBC tumor cells",
        x = "# cells",
        y = "# correctly called genotypes" #"number of calls",
        ) +
  scale_fill_manual(
    values = software_pal
  )
ggsave(
  filename="Wang2014_validated_somatic_genotype_recall.pdf",
  plot=q, 
  width=10,
  height=10
)
