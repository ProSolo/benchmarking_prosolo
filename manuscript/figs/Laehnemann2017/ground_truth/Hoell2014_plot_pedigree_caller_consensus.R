library(tidyverse)
library('RColorBrewer')
library('extrafont')
library('UpSetR')
loadfonts(device = "pdf")

sets <-   tibble(
            tool_set =
              c( "FamSeq",
                 "beagle",
                 "polymutt",
                 "FamSeq&beagle",
                 "FamSeq&polymutt",
                 "beagle&polymutt",
                 "FamSeq&beagle&polymutt"
                 )
          )

current <- tibble()

file_all <- str_c("CCS_1_W.gatk_snps.all_polymutt_beagle_FamSeqMCMC.gtisec")
all <- as_tibble( read_tsv(file_all, skip = 10, col_names = c("all"), col_types = "i") )

file_cons <- str_c("CCS_1_W.gatk_snps.cons_no_reheader_polymutt_beagle_FamSeqMCMC.gtisec")
cons <- as_tibble( read_tsv(file_cons, skip = 10, col_names = c("cons"), col_types = "i") )

sets <- sets %>% 
          add_column(all = all$all, cons = cons$cons) %>%
          mutate(diff = all - cons)

upset_cons <- fromExpression(
                sets %>% select("tool_set", "cons") %>%
                  spread(key = "tool_set", value = "cons")
                ) %>% add_column(cons = "C")

upset_diff <- fromExpression(
                sets %>% select("tool_set", "diff") %>%
                  spread(key = "tool_set", value = "diff")
                ) %>% add_column(cons = "R")

upset <- bind_rows(upset_cons, upset_diff) %>%
          mutate(cons = factor(cons),
                 beagle = as.integer(beagle),
                 FamSeq = as.integer(FamSeq),
                 polymutt = as.integer(polymutt)
                 )

cairo_pdf(file = "Hoell2014_pedigree_consensus_calling.pdf", width = 14, height = 8, onefile = F, family = "Lato")
upset(upset,
      queries = list(
                  list( query = elements,
                        params = list("cons", levels(upset$cons)[1]),
                        color = "black",
                        active = TRUE )
                ),
#      query.legend = "top",
      main.bar.color = gray(0.6),
      sets.bar.color = gray(0.6),
      mainbar.y.label = "sites with shared\ngenotype calls",
      sets.x.label = "sites with\ngenotype calls",
      set_size.show = TRUE,
      set_size.scale_max = 3000000,
      text.scale = c(4, 4, 4, 2, 4, 2),
      point.size = 3,
      mb.ratio = c(0.55, 0.45)
      )
dev.off()
