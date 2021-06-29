# packages, colors, palettes and stuff
source("../common.R")

# Read in validated somatics annotations
validated_clonals <- read_tsv(
  "nature13600-s1-table-s6-clonal.tsv"
) %>%
  mutate(
    Duplex_Freq = as.numeric( str_replace(Duplex_Freq, ",", ".") ),
    Duplex_P_val = as.numeric(
      str_replace(
        str_replace(
          str_replace(Duplex_P_val, "-[^\\d]+", "-"),
          "<", ""),
        ",", ".")
      )
    )

reclassification_motivation <- ggplot(validated_clonals) +
  geom_histogram(
    aes(
      Duplex_Freq
    ),
    bins = 50
  ) +
  scale_x_continuous(
    limits = c(0.0, 1.0),
    n.breaks = 11
  ) +
  facet_grid(
    cols = vars(class),
    rows = vars(zygosity)
  ) +
  theme_bw() +
  xlab("alternative allele frequency (duplex sequencing)") +
  ylab("clonal variant calls") 


ggsave(
  filename="Wang2014_validated_somatic_clonal_zygosity_misclassifications.pdf",
  plot=reclassification_motivation,
  width=10,
  height=5
)
