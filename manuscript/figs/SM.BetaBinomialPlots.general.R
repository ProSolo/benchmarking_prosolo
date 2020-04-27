library(VGAM)
library(ggplot2)
library(stringr)
library(RColorBrewer)

dot_to_minus <- function(x) {
  str_replace(x, '\\.', '-')
}

plot_beta_binomial <- function(cov, white, black, shape_scale) {
  p <- ggplot(data.frame(x=c(0,cov)), aes(x)) +
    stat_function(fun=dbetabinom.ab, n = cov+1, args=list(size=cov, shape1 = shape_scale * black, shape2 = shape_scale * white),
      geom="bar", fill = brewer.pal(4, "Greys")[3]) +
    theme_bw(base_size = 20, base_family = "Helvetica-Narrow") +
    labs(
      x = "alternative allele count",
      y = "probability",
      subtitle = str_c("w_t = ", white, ", b_t = ", black, ", f = ", shape_scale)
    )
  print(p)
  ggsave(str_c("Beta_binom_cov-", cov,
               "_w-", dot_to_minus(white),
               "_b-", dot_to_minus(black),
               "_f-", dot_to_minus(shape_scale),
               ".pdf"), device = cairo_pdf, width=9, height=6)
}

cov <- 20

plot_beta_binomial(cov, 2, 2, 1)
plot_beta_binomial(cov, 2, 2, .5)
plot_beta_binomial(cov, 2, 2, .25)

black_error <- 0.00000295 * cov
white <- 4
scaling <- .5

plot_beta_binomial(cov, white, black_error, scaling)


ggplot(data.frame(x=c(0,cov)), aes(x)) +
  stat_function(fun=dbetabinom.ab, n = cov+1, args=list(size=cov, shape1 = scaling * black_error, shape2 = scaling * white),
    geom="bar", fill = brewer.pal(4, "Greys")[3]) +
  theme_bw(base_size = 20, base_family = "Helvetica-Narrow") +
  coord_cartesian(ylim = c(0,.0001) ) +
  labs(
    x = "",
    y = ""
  )
ggsave(str_c("Beta_binom_cov-", cov,
             "_w-", dot_to_minus(white),
             "_b-", dot_to_minus(black_error),
             "_f-", dot_to_minus(scaling),
             "_focus_tail.pdf"), device = cairo_pdf, width=9, height=6)
