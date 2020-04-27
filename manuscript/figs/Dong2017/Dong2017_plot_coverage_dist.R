library(tidyverse)
library('RColorBrewer')
library('extrafont')
loadfonts(device = "pdf")

samples <- c("IL-11", "IL-12", "IL-1c", "Clones")

# Read in data

current <- tibble()
coverage_dists <- tibble()
for (sample in samples) {
  file <- str_c( sample, ".bps.mosdepth.global.dist.txt")
      
  current <- as_tibble( read_tsv( file,
                                  col_names = c("region", "coverage", "fraction"),
                                  col_types = "cid"
                                  ) ) %>%
              add_column( sample = sample ) %>%
              mutate( sample = factor(sample, levels = samples, ordered = TRUE) )
  
  coverage_dists <- bind_rows(coverage_dists, current)
}

linetype_values = c(3, 2, 1, 4, 5)

max_power <- 2
powers <- 10**(0:max_power)
majors <- c(1,3,5)
minors <- 1:9
major_ticks <- rep(majors, times = max_power+1)*rep(powers, each = length(majors))
minor_ticks <- rep(minors, times = max_power+1)*rep(powers, each = length(minors))

ggplot(coverage_dists %>% filter(region == "total") %>% filter(coverage != 0), aes(x = coverage, y = fraction) ) +
  coord_cartesian( ylim = c(0,1), xlim = c(1,500) ) +
  scale_x_continuous( limits = c(1,500),
                      trans = 'log10',
                      breaks = major_ticks,
                      minor_breaks = minor_ticks
                      ) +
  scale_y_continuous( limits = c(0,1),
                      breaks = seq( from = 0, to = 1, by = 1/5),
                      minor_breaks = seq( from = 0, to = 1, by = 1/10) ) +
  geom_line( aes( linetype = sample, colour = sample ) ) +
  scale_linetype_manual(
    values = linetype_values
    ) +
  scale_colour_brewer(palette = "Set1") +
  theme_bw(base_size=24, base_family="Lato")  +
  labs(x = "site coverage",
       y = "fraction of sites",
       tag = "A",
       title = "Dong 2017, whole genome")
ggsave("Dong2017_coverage_dist.pdf", device = cairo_pdf, width=10, height=6)
