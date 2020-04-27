library(VGAM)

# model for heterozygous sites
al_1 <- function(cov) {
  return(0.057378844*cov + 0.669733191)
}

al_2 <- function(cov) {
  return(0.003233912*cov + 0.399261625)
}

w <- function(cov) {
  return(0.000548761*cov + 0.540396786)
}

het_model <- function(x,cov) {
  return(w(cov)*dbetabinom.ab(x,cov,al_1(cov),al_1(cov)) + (1-w(cov))*dbetabinom.ab(x,cov,al_2(cov),al_2(cov)))
}

## generate example value vectors for het model
for (c in c(5,10,30,60)) {
  k <- 0:c
  vals <- sapply(k, het_model, cov = c)
  print(c)
  print(format(vals, nsmall = 12, scientific = FALSE))
}

# model for hom ref sites

alpha <- function(cov) {
  return(-0.000027183*cov + 0.068567471)
}

beta <- function(cov) {
  return(0.007454388*cov + 2.367486659)
}

hom_ref_model <- function(x,cov) {
  return( dbetabinom.ab(x,cov,alpha(cov),beta(cov)) )
}

## generate example value vectors for hom ref model
for (c in c(5,10,30,60)) {
  k <- 0:c
  vals <- sapply(k, hom_ref_model, cov = c)
  print(c)
  print(format(vals, nsmall = 12, scientific = FALSE))
}

# model for hom alt sites (mirror of hom ref model)
hom_alt_model <- function(x,cov) {
  return( dbetabinom.ab(x,cov,beta(cov),alpha(cov)) )
}

## generate example value vectors for hom alt model
for (c in c(5,10,30,60)) {
  k <- 0:c
  vals <- sapply(k, hom_alt_model, cov = c)
  print(c)
  print(format(vals, nsmall = 12, scientific = FALSE))
}

library(ggplot2)
library(RColorBrewer)
library(stringr)

# coverages to plot
covs <- c(20)

# plot het model

base <- ggplot(data.frame(x=c(5,100)), aes(x)) + coord_cartesian(ylim=c(0,7))
base +
  stat_function(fun=al_1, linetype=2, size = 1.2) +
  stat_function(fun=al_2, linetype=3, size = 1.2) + 
  stat_function(fun=w, linetype=4, size = 1.2)
ggsave("Lodato_het-model_parameter-scaling.pdf", device = cairo_pdf, width=9, height=6)

for (c in covs) {
  b2 <- ggplot(data.frame(x=c(0,c)), aes(x)) #+ coord_cartesian(ylim = c(0,0.1))
  b2 <- b2 + 
    stat_function(fun=het_model, n = c+1, args=list(cov=c),
                  #geom="line", colour=pal[2])
                  geom="bar", fill = brewer.pal(4, "Greys")[3]) +
    stat_function(fun=dbetabinom.ab, n = c+1, args=list(size=c,shape1=al_1(c),shape2=al_1(c)),
                  #geom="line", colour=pal[1]) + 
                  geom="line", linetype=2, size = 1.2) + 
    stat_function(fun=dbetabinom.ab, n = c+1, args=list(size=c,shape1=al_2(c),shape2=al_2(c)), 
                  #geom="line", colour=pal[3]) + 
                  geom="line", linetype=3, size = 1.2) +
    labs(tag = "B",
         title = "known heterozygous sites:\n allelic bias / dropout",
         x = "number of reads with alternative allele",
         y = "fraction of sites") +
    theme_bw(base_size = 28, base_family = "Helvetica-Narrow")
  print(b2)
  ggsave(str_c("Lodato_het-model_cov-", c, ".pdf"), device = cairo_pdf, width=11, height=7)
}

# plot hom ref model

base <- ggplot(data.frame(x=c(0,60)), aes(x)) + coord_cartesian(ylim=c(0,4.3))
base +
  stat_function(fun=alpha, linetype=2, size = 1.2) +
  stat_function(fun=beta, linetype=3, size = 1.2)
ggsave("Lodato_hom-model_parameter-scaling.pdf", device = cairo_pdf, width=9, height=6)

for (c in covs) {
  b2 <- ggplot(data.frame(x=c(0,c)), aes(x)) #+ coord_cartesian(ylim = c(0,0.1))
  b2 <- b2 + 
    stat_function(fun=hom_ref_model, n = c+1, args=list(cov=c), geom="bar", fill = brewer.pal(4, "Greys")[3]) + 
    labs(tag = "A",
         title = "known homozygous reference sites:\n alternative allele coverage",
         x = "number of reads with alternative allele",
         y = "fraction of sites") +
    theme_bw(base_size = 28, base_family = "Helvetica-Narrow")
  print(b2)
  ggsave(str_c("Lodato_hom-ref-model_cov-", c, ".pdf"), device = cairo_pdf, width=11, height=7)
}

# plot hom alt model

for (c in covs) {
  b2 <- ggplot(data.frame(x=c(0,c)), aes(x)) #+ coord_cartesian(ylim = c(0,0.1))
  b2 <- b2 + 
    stat_function(fun=hom_alt_model, n = c+1, args=list(cov=c), geom="bar", fill = brewer.pal(4, "Greys")[3]) +
    labs(tag = "C",
         title = "known hom alt sites: ref coverage",
         x = "reads with alternative nucleotide",
         y = "fraction of sites") +
    theme_bw(base_size = 28, base_family = "Helvetica-Narrow")
  print(b2)
  ggsave(str_c("Lodato_hom-alt-model_cov-", c, ".pdf"), device = cairo_pdf, width=11, height=6)
}

