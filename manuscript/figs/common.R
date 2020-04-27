library(tidyverse)
library('RColorBrewer')
library('ggthemes')
library('extrafont')
loadfonts(device = "pdf")

linetype_pal = c(
  "MonoVar default" = 1,
  "MonoVar no consensus" = 3,
  "ProSolo default" = 1,
  "ProSolo imputation" = 3,
  "SCAN-SNV sensitive" = 3,
  "SCcaller default dbsnp" = 1,
  "SCcaller sensitive dbsnp" = 3,
  "SCIPhI default" = 1,
  "SCIPhI sensitive" = 3
  )

# Based on Okabe and Ito colorblind-safe palette (https://jfly.uni-koeln.de/color/)
#' Tint a color, i.e. make it lighter by a factor
#' 
#' @description Tint a given color by a given factor
#' 
#' @param col HEX valued color to lighten up.
#' @param factor Factor to use for tinting. A higher factor will produce a lighter color.
tint <- function(col, factor) {
  rgb_col <- col2rgb(col)
  rgb(t( ( (255 - rgb_col ) * factor ) + rgb_col), maxColorValue=255)
}

factor <- 1/3

Black <- colorblind_pal()(8)[[1]]
light_Black <- tint(Black, factor)

Vermillion_Red <- colorblind_pal()(8)[[7]]
light_Vermillion_Red <- tint(Vermillion_Red, factor)

Blue <- colorblind_pal()(8)[[6]]
light_Blue <- tint(Blue, factor)

reddish_Purple <- colorblind_pal()(8)[[8]]
light_reddish_Purple <- tint(reddish_Purple, factor)

Yellow <- colorblind_pal()(8)[[5]]

# unused Okabe Ito colors
# Sky_Blue <- colorblind_pal()(8)[[3]]
# Orange <- colorblind_pal()(8)[[2]]
# bluish_Green <- colorblind_pal()(8)[[4]]

# a look at the selected colors
# n_col <- 9
# dummy_counts <- tibble(color = factor(seq_len(n_col)), count = rep(4, times=n_col))
# 
# bp <- ggplot(dummy_counts, aes(color, count)) +
#         geom_col(aes(fill=color) ) +
#         scale_fill_manual(
#           values=c( Black, light_Black,
#                 Vermillion_Red, light_Vermillion_Red,
#                 Blue, light_Blue,
#                 reddish_Purple, light_reddish_Purple,
#                 Yellow 
#                 )
#          )
# bp

# requires install of colorblindr package, not on conda-forge, but displays
# how people with different colorblindness perceive the colors
# library(colorblindr)
# cvd_grid(bp)

software_pal <- c(
  "MonoVar default" = Black,
  "MonoVar no consensus" = light_Black,
  "ProSolo default" = Vermillion_Red,
  "ProSolo imputation" = light_Vermillion_Red,
  "SCAN-SNV sensitive" = Yellow,
  "SCcaller default dbsnp" = Blue,
  "SCcaller sensitive dbsnp" = light_Blue,
  "SCIPhI default" = reddish_Purple,
  "SCIPhI sensitive" = light_reddish_Purple
)

shape_pal <-
  c("1e-06" = 17,
    "1e-05" = 4,
    "1e-04" = 8,
    "0.001" = 1,
    "0.005" = 7,
    "0.01" = 2,
    "0.02" = 9,
    "0.05" = 13,
    "0.1" = 3,
    "0.2" = 10,
    "1" = 20
  )
