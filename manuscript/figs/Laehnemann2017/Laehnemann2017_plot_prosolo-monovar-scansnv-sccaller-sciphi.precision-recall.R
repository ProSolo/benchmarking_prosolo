# packages, colors, palettes and stuff
source("../common.R")

cells <- c("PAG1", "PAG2", "PAG5", "PAG9", "PAG10")

# Read in MonoVar stats

# TODO: add "0-1" once computed
ps <- c( "0-000001", "0-00001", "0-0001", "0-001", "0-01", "0-02", "0-05", "0-1", "0-2", "1")
cs <- c("0", "1") # indicator whether consensus filtering was on in MonoVar

current <- tibble()
positives_negatives <- tibble()
for (cell in cells) {
  for ( c in cs ) {
    for ( p in ps ) {
      file <- str_c("monovar/", cell, ".single_cells.", p, ".c", c, ".monovar.positives_negatives.tsv")
      
      if ( c == 1 ) {
        mode = "MonoVar default"
      } else {
        mode  = "MonoVar no consensus"
      }
      
      current <- as_tibble( read_tsv(file) ) %>%
                  add_column( software = mode,
                              cell = cell,
                              filter = if_else(p == "1", 1.0, as.numeric( str_replace(p, '-', '.') )) )
      
      positives_negatives <- bind_rows(positives_negatives, current)
    }
  }
}

# Read in SCcaller stats

alphas <- c("0-01", "0-05", "1-0" ) # different alpha values for filtering SCcaller results, 1.0 means no filtering
dbsnp <- c("", ".dbsnp") # against custom het sites or known dbsnp sites
low_cov <- c("", ".low-cov") # with repo settings or low-cov settings similar to publication

for (cell in cells) {
  for ( d in dbsnp ) {
    for ( m in low_cov ) {
      for ( a in alphas ) {
          file <- str_c("sccaller/", cell, d, m, ".sccaller.varcall.cutoff.fil-amp-err_alpha_", a, ".positives_negatives.tsv")
          
          current <- as_tibble( read_tsv(file) ) %>%
                      add_column( software = case_when(
                                      d != "" & m == "" ~ "SCcaller default dbsnp",
                                      d != "" & m != "" ~ "SCcaller sensitive dbsnp",
                                      d == "" & m == "" ~ "SCcaller default bulk",
                                      d == "" & m != "" ~ "SCcaller sensitive bulk"),
                                  cell = cell,
                                  filter = if_else(a == "1-0", 1.0, as.numeric( str_replace(a, '-', '.') ) ) )
          
          positives_negatives <- bind_rows(positives_negatives, current)
      }
    }
  }
}

# Read in ProSolo stats

fdrs <- c( "0-000001", "0-00001", "0-0001", "0-001", "0-005", "0-01", "0-02", "0-05", "0-1", "0-2", "1")
modes <- c( "", ".min_sc_cov_1" )

for (cell in cells) {
  for ( f in fdrs ) {
    for ( m in modes ) {
      file <- str_c( "prosolo/fdr_alt-presence/", cell, ".PNG.alt_sites_only.fdr_", f, "_alt-presence.prosolo", m , ".positives_negatives.ML_set.tsv" )
      
      if ( m == ".min_sc_cov_1" ) {
        mode = "ProSolo default"
      } else {
        mode = "ProSolo imputation"
      }
    
      current <- as_tibble( read_tsv(file) ) %>%
                  add_column( software = mode,
                              cell = cell,
                              filter = if_else(f == "1", 1.0, as.numeric( str_replace(f, '-', '.') ) ) )
      
      positives_negatives <- bind_rows(positives_negatives, current)
    }
  }
}

# Read in SCAN-SNV stats

fdrs <- c( "0-000001", "0-00001", "0-0001", "0-001", "0-005", "0-01", "0-02", "0-05", "0-1", "0-2", "1-0")

for (cell in cells) {
  for ( f in fdrs ) {
    file <- str_c( "scansnv/", cell, ".alt.", f, ".positives_negatives.tsv" )
    
    mode = "SCAN-SNV sensitive"
    
    current <- as_tibble( read_tsv(file) ) %>%
                add_column( software = mode,
                            cell = cell,
                            filter = as.numeric( str_replace(f, '-', '.') )
                            )
    
    positives_negatives <- bind_rows(positives_negatives, current)
  }
}

# Read in SciPhi stats

modes <- c("default", "sensitive")

for (cell in cells) {
  for (m in modes) {
    if (m == "default") {
      file <- str_c("sciphi/", cell, ".single_cells.sciphi.positives_negatives.tsv")
      mode = "SCIPhI default"
    } else if (m == "sensitive") {
      file <- str_c("sciphi/", cell, ".single_cells.400000_iterations.sensitive.sciphi.positives_negatives.tsv")
      mode = "SCIPhI sensitive"
    }
    
    current <- as_tibble( read_tsv(file) ) %>%
                add_column( software = mode,
                            cell = cell,
                            filter = 1.0)
    
    positives_negatives <- bind_rows(positives_negatives, current)
  }
}

metrics <- positives_negatives %>%
              mutate( recall = TP / ( TP + FN + P_not_called),
                      FPR = FP / (FP + TN + N_not_called),
                      FDR = FP / (FP + TP),
                      precision = TP / ( TP + FP ),
                      F1 = 2 * ( (precision * recall) / (precision + recall) ),
                      F0_5 = 1.25 * ( (precision * recall ) / ( 0.25 * precision + recall ) ),
                      filter_factor = factor(filter)
                    )

avg_per_software_and_parameter <- metrics %>%
                  complete( software, filter ) %>%
                  mutate(filter = factor(filter, ordered = TRUE)) %>%
                  group_by(software, filter) %>%
                  summarise( avg_prec = mean(precision),
                             min_prec = min(precision),
                             max_prec = max(precision),
                             avg_FPR = mean(FPR),
                             min_FPR = min(FPR),
                             max_FPR = max(FPR),
                             avg_FDR = mean(FDR),
                             min_FDR = min(FDR),
                             max_FDR = max(FDR),
                             avg_rec = mean(recall),
                             min_rec = min(recall),
                             max_rec = max(recall) )

max_rec_for_prec_above_99 <- avg_per_software_and_parameter %>%
    filter(avg_prec > 0.99) %>%
    select(software, filter, avg_prec, avg_rec) %>%
    separate(software, sep = " ", into = c("software", "mode", "mode2")) %>%
    filter(mode != "imputation") %>%
    group_by(software) %>%
    summarise(max(avg_rec))
write_tsv(max_rec_for_prec_above_99, "Laehnemann2017_max_avg_rec_for_prec_above_99.tsv")

fils = levels(avg_per_software_and_parameter$filter)

axis_limits <- c(0.0000001,1)

ggplot(avg_per_software_and_parameter %>%
         filter(software == "ProSolo default" | software == "ProSolo imputation") %>%
         mutate(num_filter = as.numeric(levels(filter)[filter]) )
  ) + 
  geom_point(aes(x = num_filter, y = avg_FDR, color = software, shape = filter)) +
  geom_abline(intercept = 0) +
  scale_x_continuous(
    trans = "log10",
    limits=axis_limits
  ) +
  scale_y_continuous(
    trans = "log10",
    limits=axis_limits
  ) +
  scale_color_manual(
    values = software_pal
  ) +
  scale_shape_manual(
    name = "theoretical FDR",
    values = shape_pal
  ) +
  theme_bw(base_size=24, base_family="Lato")  +
  labs(
    x = "theoretical FDR",
    y = "ground truth FDR",
    tag = "B",
    title = "5 granulocytes"
  )
ggsave("Laehnemann2017_prosolo__FDR_ground_truth_vs_theoretical.pdf", device = cairo_pdf, width=11, height=7.5)
 

ggplot(metrics %>%
        filter(software != "SCcaller default bulk", software != "SCcaller sensitive bulk"),
  aes(x = filter, y = F1) ) +
  coord_cartesian( ylim = c(0,1) ) +
  geom_line( aes(group = interaction(cell, software), color = software, linetype = software)) +
  geom_point( aes(shape = filter_factor, color = software, fill = software), size = 4) +
  scale_x_continuous(
    trans = 'log10'
  ) +
  scale_shape_manual(
    guide = "none",
    values = shape_pal
    ) +
  scale_colour_manual(
    guide = "none",
    values = software_pal
    ) +
  scale_fill_manual(
    guide = "none",
    values = software_pal
    ) +
  scale_linetype_manual(
    guide = "none",
    values = linetype_pal
    ) +
#  guides(
#    shape = guide_legend( title = "threshold used") ) +
  theme_bw(base_size=24, base_family="Lato")  +
  theme( axis.text.x = element_text( angle=45, hjust = 1) ) +
  facet_wrap(~cell, nrow=3) +
  labs(x = "threshold used",
       y = "F_1 score",
       tag = "B",
       title = "5 granulocytes")
ggsave("Laehnemann2017_prosolo-monovar-scansnv-sccaller-sciphi_F1_plot.pdf", device = cairo_pdf, width=14, height=20)

ggplot(avg_per_software_and_parameter %>%
    filter(software != "SCcaller default bulk", software != "SCcaller sensitive bulk"),
    aes(x = avg_rec, y = avg_prec) ) +
  coord_cartesian( ylim = c(0,1), xlim = c(0,1) ) +
  geom_line( aes( color = software, linetype = software ) ) +
  scale_linetype_manual(
    values = linetype_pal,
    labels = software_labels
    ) +
  geom_point( aes(shape = filter, color = software, fill = software), size = 4) +
  scale_colour_manual(
    values = software_pal,
    labels = software_labels
    ) +
  scale_fill_manual(
    values = software_pal,
    labels = software_labels
    ) +
  scale_shape_manual(
    values = shape_pal
    ) +
  guides(
    linetype = guide_legend( keywidth = 3, keyheight = 1.5),
    shape = guide_legend( title = "threshold used")
    ) +
  theme_bw(base_size=24, base_family="Lato")  +
  labs(x = "average recall",
       y = "average precision",
       tag = "B",
       title = "5 granulocytes")
ggsave("Laehnemann2017_prosolo-monovar-scansnv-sccaller-sciphi_precision-recall-plot.pdf", device = cairo_pdf, width=10, height=7.5)

# focus on main tool area
ggplot(avg_per_software_and_parameter %>%
    filter(software != "SCcaller default bulk", software != "SCcaller sensitive bulk"),
    aes(x = avg_rec, y = avg_prec) ) +
  coord_cartesian( ylim = c(0.875,1), xlim = c(0,0.35) ) +
  geom_line( aes( color = software, linetype = software ) ) +
  scale_linetype_manual(
    values = linetype_pal,
    labels = software_labels,
    guide = "legend"
    ) +
  geom_point( aes(shape = filter, color = software, fill = software), size = 4) +
  scale_shape_manual(
    guide = "legend",
    values = shape_pal
    ) +
  scale_colour_manual(
    guide = "legend",
    labels = software_labels,
    values = software_pal
    ) +
  scale_fill_manual(
    guide = "legend",
    labels = software_labels,
    values = software_pal
    ) +
  guides( 
    linetype = guide_legend( keywidth = 3, keyheight = 1.5),
    shape = guide_legend( title = "threshold used")
    ) +
  theme_bw(base_size=24, base_family="Lato") +
  labs(x = "average recall",
       y = "average precision",
       tag = "B",
       title = "5 granulocytes")
ggsave("Laehnemann2017_prosolo-monovar-scansnv-sccaller-sciphi_precision-recall-plot_focus-top-left.pdf", device = cairo_pdf, width=10, height=7.5)
