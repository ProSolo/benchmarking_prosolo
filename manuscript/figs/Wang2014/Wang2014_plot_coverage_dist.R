library(tidyverse)
library('RColorBrewer')
library('extrafont')
loadfonts(device = "pdf")

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
    "TNBC-n9",
    "TNBC-Pop-Tumor",
    "TNBC-Pop-Normal"
  ),
  batch = c(
    rep("tumor single cell", times = 16),
    rep("normal single cell", times = 16),
    rep("bulk", times = 2)
  )
)

# Read in data

current <- tibble()
coverage_dists <- tibble()
for (sample in cells$name) {
  file <- str_c( sample, ".bps.TruSeq_exome_targeted_regions.mosdepth.region.dist.txt")
      
  batch = cells %>% filter(name == sample) %>% pull(batch)
  current <- as_tibble(
      read_tsv(
        file,
        col_names = c("region", "coverage", "fraction"),
        col_types = "cid"
      )
    ) %>%
    add_column(
      sample = sample,
      batch = batch
    ) %>%
    mutate(
      sample = factor(sample, levels = cells$name, ordered = TRUE)
    )
  
  coverage_dists <- bind_rows(coverage_dists, current)
}

max_power <- 2
powers <- 10**(0:max_power)
majors <- c(1,3,5)
minors <- 1:9
major_ticks <- rep(majors, times = max_power+1)*rep(powers, each = length(majors))
minor_ticks <- rep(minors, times = max_power+1)*rep(powers, each = length(minors))

# we need 32 somewhat distinguishable colors, so we go for four 9-value
# sequential palettes without their lightest color -- palettes are only
# colorblind-safe within themselves, but together with linetypes this
# hopefully works for everybody
palette <- c(
  brewer.pal(9, "YlGn")[2:9],
  brewer.pal(9, "YlOrRd")[2:9],
  brewer.pal(9, "YlGnBu")[2:9],
  brewer.pal(9, "RdPu")[2:9],
  brewer.pal(4, "Greys")[3:4]
)


ggplot(coverage_dists %>% filter(region == "total") %>% filter(coverage != 0), aes(x = coverage, y = fraction) ) +
  coord_cartesian( ylim = c(0,1), xlim = c(1,200) ) +
  scale_x_continuous( limits = c(1,200),
                      trans = 'log10',
                      breaks = major_ticks,
                      minor_breaks = minor_ticks
                      ) +
  scale_y_continuous( limits = c(0,1),
                      breaks = seq( from = 0, to = 1, by = 1/5),
                      minor_breaks = seq( from = 0, to = 1, by = 1/10) ) +
  geom_line( aes( linetype = batch, colour = sample ) ) +
  scale_color_manual(
    values = palette,
    guide = guide_legend(
      ncol = 4
    )
  ) +
  theme_bw(
    base_size=24,
    base_family="Lato"
  )  +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical"
  ) +
  labs(x = "site coverage",
       y = "fraction of sites",
       tag = "C",
       title = "Wang 2014, whole exome")
ggsave("Wang2014_coverage_dist.pdf", device = cairo_pdf, width=12, height=11)
