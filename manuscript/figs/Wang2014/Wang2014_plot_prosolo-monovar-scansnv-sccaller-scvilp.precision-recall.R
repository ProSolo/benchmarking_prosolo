# packages, colors, palettes and stuff
source("../common.R")

cells <- tibble(
  name = c(
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
    "h8",
    "TNBC-n10",
    "TNBC-n11",
    "TNBC-n12",
    "TNBC-n13",
    "TNBC-n14",
    "TNBC-n15",
    "TNBC-n16",
    "TNBC-n1",
    "TNBC-n2",
    "TNBC-n3",
    "TNBC-n4",
    "TNBC-n5",
    "TNBC-n6",
    "TNBC-n7",
    "TNBC-n8",
    "TNBC-n9"
  ),
  batch = c(
    rep("tumor", times = 16),
    rep("normal", times = 16)
  )
)

get_cells_batch <- function(cells, cell) {
  cells %>% filter(name == cell) %>% pull(batch)
}

get_batch_ground_truth <- function(batch) {
  if (batch == "normal") {
    ground_truth = ".tumor_ground_truth."
  } else if (batch == "tumor") {
    ground_truth = ".normal_ground_truth."
  } else {
    stop("Unknown batch type.")
  }
  ground_truth
}

# Read in MonoVar stats

# TODO: add "0-1" once computed
ps <- c( "0-000001", "0-00001", "0-0001", "0-001", "0-01", "0-02", "0-05", "0-1", "0-2", "1")
cs <- c("0", "1") # indicator whether consensus filtering was on in MonoVar

current <- tibble()
positives_negatives <- tibble()
for (cell in cells[["name"]]) {
  for ( c in cs ) {
    for ( p in ps ) {
      batch <- get_cells_batch(cells, cell)
      ground_truth = get_batch_ground_truth(batch)
      file <- str_c("monovar/", cell, ".", batch, "_cells", ground_truth, p, ".c", c, ".monovar.alt_calls.awk_positives_negatives.tsv")
      
      if ( c == 1 ) {
        mode = "MonoVar default"
      } else {
        mode  = "MonoVar no consensus"
      }
       
      current <- as_tibble( read_tsv(file) ) %>%
                  add_column( software = mode,
                              cell = cell,
                              batch = batch,
                              filter = if_else(p == "1", 1.0, as.numeric( str_replace(p, '-', '.') )) )
      
      positives_negatives <- bind_rows(positives_negatives, current)
    }
  }
}

# Read in SCcaller stats

alphas <- c("0-01", "0-05", "1-0" ) # different alpha values for filtering SCcaller results, 1.0 means no filtering
dbsnp <- c("dbsnp") # against custom het sites or known dbsnp sites
low_cov <- c("", ".low-cov") # with repo settings or low-cov settings similar to publication

for (cell in cells[["name"]]) {
  for ( d in dbsnp ) {
    for ( m in low_cov ) {
      for ( a in alphas ) {
          batch <- get_cells_batch(cells, cell)
          ground_truth = get_batch_ground_truth(batch)
          file <- str_c("sccaller/", cell, ground_truth, d, m, ".sccaller.varcall.cutoff.fil-amp-err_alpha_", a, ".alt_calls.awk_positives_negatives.tsv")
          
          current <- as_tibble( read_tsv(file) ) %>%
                      add_column( software = case_when(
                                      d != "" & m == "" ~ "SCcaller default dbsnp",
                                      d != "" & m != "" ~ "SCcaller sensitive dbsnp",
                                      d == "" & m == "" ~ "SCcaller default bulk",
                                      d == "" & m != "" ~ "SCcaller sensitive bulk"),
                                  cell = cell,
                                  batch = batch,
                                  filter = if_else(a == "1-0", 1.0, as.numeric( str_replace(a, '-', '.') ) ) )
          
          positives_negatives <- bind_rows(positives_negatives, current)
      }
    }
  }
}

# Read in ProSolo stats

fdrs <- c( "0-000001", "0-00001", "0-0001", "0-001", "0-005", "0-01", "0-02", "0-05", "0-1", "0-2", "0-5", "1")
modes <- c( ".min_sc_cov_1" )

for (cell in cells[["name"]]) {
  for ( f in fdrs ) {
    for ( m in modes ) {
      batch <- get_cells_batch(cells, cell)
      ground_truth = get_batch_ground_truth(batch)
      file <- str_c( "prosolo/fdr_alt-presence/", cell, ".TNBC-Pop-", str_to_sentence(batch), ground_truth, "alt_sites_only.fdr_", f, "_alt-presence.prosolo", m , ".alt_calls.awk_positives_negatives.ML_set.tsv" )
      
      if ( m == ".min_sc_cov_1" ) {
        mode = "ProSolo default"
      } else {
        mode = "ProSolo imputation"
      }
      
      current <- as_tibble( read_tsv(file) ) %>%
                  add_column( software = mode,
                              cell = cell,
                              batch = batch,
                              filter = if_else(f == "1", 1.0, as.numeric( str_replace(f, '-', '.') ) )
                              )
      
      positives_negatives <- bind_rows(positives_negatives, current)
    }
  }
}

## Read in SCAN-SNV stats
#
#fdrs <- c( "0-000001", "0-00001", "0-0001", "0-001", "0-005", "0-01", "0-02", "0-05", "0-1", "0-2", "1-0")
#
#for (cell in cells) {
#  for ( f in fdrs ) {
#    file <- str_c( "scansnv/", cell, ".alt.", f, ".positives_negatives.tsv" )
#    
#    mode = "SCAN-SNV sensitive"
#    
#    current <- as_tibble( read_tsv(file) ) %>%
#                add_column( software = mode,
#                            cell = cell,
#                            filter = as.numeric( str_replace(f, '-', '.') )
#                            )
#    
#    positives_negatives <- bind_rows(positives_negatives, current)
#  }
#}
#
# Read in SciPhi stats

modes <- c("default")
fdrs <- c( "0-000001", "0-00001", "0-0001", "0-001", "0-005", "0-01", "0-02", "0-05", "0-1", "0-2", "0-5", "1")

for (cell in cells[["name"]]) {
  for (m in modes) {
    for (f in fdrs) {
      batch <- get_cells_batch(cells, cell)
      ground_truth = get_batch_ground_truth(batch)
      file <- str_c("sciphi/fdr_alt-presence/", cell, ".", batch, "_cells", ground_truth, "80000-40000_iterations.default.sciphi.alt_prob.fdr_", f, "_alt-presence.alt_calls.awk_positives_negatives.ML_set.tsv")
      mode = "SCIPhI default"
      
      current <- as_tibble( read_tsv(file) ) %>%
                  add_column( software = mode,
                              cell = cell,
                              batch = batch,
                              filter = if_else(f == "1", 1.0, as.numeric( str_replace(f, '-', '.') ) )
                              )
      
      positives_negatives <- bind_rows(positives_negatives, current)
    }
  }
}

metrics <- positives_negatives %>%
              mutate( recall = TP / ( TP + FN + P_not_called),
                      FPR = FP / (FP + TN + N_not_called),
                      FDR = FP / (FP + TP),
                      precision = TP / ( TP + FP ),
                      F1 = 2 * ( (precision * recall) / (precision + recall) ),
                      F0_5 = 1.25 * ( (precision * recall ) / ( 0.25 * precision + recall ) ),
                      filter =
                        factor(filter,
                              levels = c( 0.000001, 0.00001, 0.0001, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
                        )
                    )

avg_per_software_batch_and_parameter <- metrics %>%
                  complete( software, filter, batch ) %>%
                  mutate(filter = factor(filter, ordered = TRUE)) %>%
                  group_by(software, filter, batch) %>%
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

for (b in unique( cells %>% pull(batch) ) ) {
  max_rec_for_prec_above_99 <- avg_per_software_batch_and_parameter %>%
    filter(batch == b) %>%
    filter(avg_prec > 0.99) %>%
    select(software, filter, avg_prec, avg_rec) %>%
    separate(software, sep = " ", into = c("software", "mode", "mode2")) %>%
    group_by(software) %>%
    summarise(max(avg_rec))
  write_tsv(max_rec_for_prec_above_99, str_c("Wang2014_max_avg_rec_for_prec_above_99_", b, ".tsv") )
}

fils = levels(avg_per_software_batch_and_parameter$filter)

axis_limits <- c(0.0000001,1)
reduced_software_pal <- software_pal %>%
  enframe() %>% 
  filter(
    name == "ProSolo default" |
      name == "ProSolo imputation" |
      name == "SCIPhI default" |
      name == "SCIPhI sensitive"
    ) %>%
  deframe()

ggplot(avg_per_software_batch_and_parameter %>%
         filter(
           software == "ProSolo default" | 
           software == "SCIPhI default"
         ) %>%
         mutate(
           num_filter = as.numeric(levels(filter)[filter]),
           software =
             factor(software,
                    levels = c("ProSolo default", "ProSolo imputation", "SCIPhI default", "SCIPhI sensitive")
             )
         )
  ) + 
  geom_point(aes(x = num_filter, y = avg_FDR, color = software, shape = filter), size = 4) +
  geom_segment(aes(x = 0.000001, y = 0.000001, xend = 1, yend = 1)) +
  scale_x_continuous(
    trans = "log10",
    limits=axis_limits
  ) +
  scale_y_continuous(
    trans = "log10",
    limits=axis_limits
  ) +
  scale_color_manual(
    values = reduced_software_pal,
    labels = software_labels,
    drop = TRUE
  ) +
  scale_shape_manual(
    name = "theoretical FDR",
    values = shape_pal,
    drop = FALSE
  ) +
  facet_grid(
    cols = vars(batch)
  ) +
  theme_bw(base_size=24, base_family="Lato")  +
  labs(
    x = "theoretical FDR",
    y = "ground truth FDR",
    tag = "c",
    title = "TNBC"
  )
ggsave(
  "Wang2014_prosolo-sciphi_FDR_ground_truth_vs_theoretical.pdf",
  device = cairo_pdf,
  width=16,
  height=7.5
)

basic_software_pal <- software_pal %>%
  enframe() %>%
  filter(
    name %in% basic_software_levels
  ) %>%
  deframe()

basic_linetype_pal <- linetype_pal %>%
  enframe() %>%
  filter(
    name %in% basic_software_levels
  ) %>%
  deframe()

ggplot(avg_per_software_batch_and_parameter %>%
    filter(
      software != "SCcaller default bulk",
      software != "SCcaller sensitive bulk"
      ) %>%
    mutate(software = factor(software, levels = basic_software_levels)),
    aes(x = avg_rec, y = avg_prec) ) +
  coord_cartesian( ylim = c(0,1), xlim = c(0,1) ) +
  geom_line( aes( color = software, linetype = software ) ) +
  scale_linetype_manual(
    values = basic_linetype_pal,
    labels = software_labels,
    drop = FALSE
    ) +
  geom_point( aes(shape = filter, color = software, fill = software), size = 4) +
  scale_colour_manual(
    values = basic_software_pal,
    labels = software_labels,
    drop = FALSE
    ) +
  scale_fill_manual(
    values = basic_software_pal,
    labels = software_labels,
    drop = FALSE
    ) +
  scale_shape_manual(
    values = shape_pal,
    drop = FALSE
    ) +
  guides(
    linetype = guide_legend( keywidth = 3, keyheight = 1.5),
    shape = guide_legend( title = "threshold used")
    ) +
  facet_grid(
    cols = vars(batch)
  ) +
  theme_bw(base_size=24, base_family="Lato")  +
  labs(x = "average recall",
       y = "average precision",
       tag = "c",
       title = "TNBC")
ggsave(
  "Wang2014_prosolo-monovar-scansnv-sccaller-scvilp_precision-recall-plot.pdf",
  device = cairo_pdf,
  width=16,
  height=8.5
)

# focus on main tool area
ggplot(avg_per_software_batch_and_parameter %>%
    filter(software != "SCcaller default bulk", software != "SCcaller sensitive bulk") %>%
    mutate(software = factor(software, levels = basic_software_levels)),
    aes(x = avg_rec, y = avg_prec) ) +
  coord_cartesian( ylim = c(0.9,1), xlim = c(0.1,0.7) ) +
  geom_line( aes( color = software, linetype = software ) ) +
#  geom_errorbar( aes( ymin = min_prec, ymax = max_prec, color = software) ) +
#  geom_errorbar( aes( xmin = min_rec, xmax = max_rec, color = software ) ) +
  scale_linetype_manual(
    values = basic_linetype_pal,
    labels = software_labels,
    drop = FALSE
    ) +
  geom_point( aes(shape = filter, color = software, fill = software), size = 4) +
  scale_colour_manual(
    values = basic_software_pal,
    labels = software_labels,
    drop = FALSE
    ) +
  scale_fill_manual(
    values = basic_software_pal,
    labels = software_labels,
    drop = FALSE
    ) +
  scale_shape_manual(
    values = shape_pal,
    drop = FALSE
    ) +
  guides(
    linetype = guide_legend( keywidth = 3, keyheight = 1.5),
    shape = guide_legend( title = "threshold used")
    ) +
  facet_grid(
    cols = vars(batch)
  ) +
  theme_bw(base_size=24, base_family="Lato") +
  labs(x = "average recall",
       y = "average precision",
       tag = "c",
       title = "TNBC")
ggsave(
  "Wang2014_prosolo-monovar-scansnv-sccaller-scvilp_precision-recall-plot_focus-tools.pdf",
  device = cairo_pdf,
  width=16,
  height=8.5
)

for (b in unique( cells %>% pull(batch) ) ) {
  if (b == "tumor") {
    ylim <- c(0.91,1)
    xlim <- c(0.1,0.41)
  } else if (b == "normal") {
    ylim <- c(0.95,1)
    xlim <- c(0.1,0.7)
  }
  p <- ggplot(
      metrics %>%
      filter(batch == b) %>%
      mutate(software = factor(software, levels = basic_software_levels)),
      aes(x = recall, y = precision) ) +
    coord_cartesian( ylim = ylim, xlim = xlim ) +
    geom_line( aes( color = software, linetype = software) ) +
    scale_linetype_manual(
      values = basic_linetype_pal,
      labels = software_labels,
      drop = FALSE
      ) +
    geom_point( aes(shape = filter, color = software, fill = software), size = 4) +
    scale_colour_manual(
      values = basic_software_pal,
      labels = software_labels,
      drop = FALSE
      ) +
    scale_fill_manual(
      values = basic_software_pal,
      labels = software_labels,
      drop = FALSE
      ) +
    scale_shape_manual(
      values = shape_pal,
      drop = FALSE
      ) +
    guides(
      linetype = guide_legend( keywidth = 3, keyheight = 1.5),
      shape = guide_legend( title = "threshold used")
      ) +
    theme_bw(base_size=24, base_family="Lato") +
    facet_wrap(
      facets=vars(cell),
      nrow = 4
    ) +
    labs(x = "recall",
         y = "precision"
    )
  p
  ggsave(
    str_c("Wang2014_prosolo-monovar-scansnv-sccaller-scvilp_per-cell_precision-recall-plot_focus-tools_", b, ".pdf"),
    plot = p,
    device = cairo_pdf,
    width = 14,
    height = 10
  )
}
