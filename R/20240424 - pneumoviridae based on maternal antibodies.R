setwd("~/R/Larman_Lab")
library(tidyverse)
library(RColorBrewer)
library(gtable)
library(grid)
library(ggplot2)
library(ggbeeswarm)
source("ArboScan/common_functions.R")
library(data.table)

########
# Data Prep
#######
library(readxl)
sample_meta <- read_excel("BBC/sample_meta.xls")

sample_meta <- sample_meta %>% mutate(age = difftime(time1 = Plasma_Date, time2 = DOD, units = "days")) %>%
  mutate(years = (age %>% as.numeric())/365.25) %>% mutate(Plasma_type = ifelse(Plasma_type=="Postnatal","child",Plasma_type))

BBC_VARscores <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240424 - BBC VARscores - maternal abs removed.rds")

BBC_VARscores <- BBC_VARscores %>% dplyr::rename(vir_score = value) %>%
  left_join(sample_meta, multiple = "all")


hsv_VAR <- BBC_VARscores %>%
  filter(taxon_genus %in% c("Orthopneumovirus", "Metapneumovirus")) %>%
  filter(taxon_species != "Human respiratory syncytial virus") %>%
  mutate(taxon_species = as.factor(taxon_species)) %>%
  arrange(taxon_species)


Maternal_pos <- hsv_VAR %>% ungroup() %>%
  pivot_wider(id_cols = c(ID, Plasma_type), names_from = taxon_species,
              values_from = vir_score, values_fn = max) %>%
  filter(Plasma_type == "cord")

Maternal_pos$`Human respiratory syncytial virus A` %>% quantile()

Maternal_pos <- Maternal_pos %>%
  mutate(`Mom RSVA quartile` = case_when(
    `Human respiratory syncytial virus A` < 7.413402 ~ "1st quartile",
    `Human respiratory syncytial virus A` < 8.975688 ~ "2nd quartile",
    `Human respiratory syncytial virus A` < 11.119927 ~ "3rd quartile",
    T ~ "4th quartile",
    
  )) %>%
  dplyr::select(`ID`, `Mom RSVA quartile`)


dotplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
               left_join(Maternal_pos) %>% filter(taxon_species == "Human respiratory syncytial virus A") %>%
               filter(!is.na(`Mom RSVA quartile`)),
             aes(x = `Mom RSVA quartile`, y = vir_score,
                 color = `Mom RSVA quartile`),
             position = position_quasirandom(),
             size = .1,
             alpha = .5) +
  geom_hline(yintercept = 1.936591) +
  ylab("Child RSV A VARscore") +
  theme_light()+
  scale_color_manual(values = c("black", "#005A32", "#41AB5D", "#A1D99B")) +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())

wilcox.test((hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
              left_join(Maternal_pos) %>% filter(taxon_species == "Human respiratory syncytial virus A") %>%
              filter(!is.na(`Mom RSVA quartile`)) %>% filter(`Mom RSVA quartile` == "1st quartile"))$vir_score,
            (hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
               left_join(Maternal_pos) %>% filter(taxon_species == "Human respiratory syncytial virus A") %>%
               filter(!is.na(`Mom RSVA quartile`)) %>% filter(`Mom RSVA quartile` == "4th quartile"))$vir_score) # not significantly different

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - RSV A VARscore dotplot by maternal RSV A antibodies.jpeg",
       dotplot,
       scale = 1,
       height = 2,
       width = 3.5,
       units = "in",
       dpi = 300)

positive_rate <- hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
  left_join(Maternal_pos) %>% filter(taxon_species == "Human respiratory syncytial virus A") %>%
  filter(!is.na(`Mom RSVA quartile`)) %>%
  mutate(positive = vir_score > 1.936591) %>%
  group_by(`Mom RSVA quartile`) %>%
  summarise(n_positive = sum(positive), n_total = n()) %>%
  mutate(proportion = n_positive / n_total)

fisher.test(matrix(data = c(115, 127,
                            136, 106), nrow = 2))



dotplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
               left_join(Maternal_pos) %>% filter(taxon_species == "Human respiratory syncytial virus B") %>%
               filter(!is.na(`Mom RSVA quartile`)),
             aes(x = `Mom RSVA quartile`, y = vir_score,
                 color = `Mom RSVA quartile`),
             position = position_quasirandom(),
             size = .1,
             alpha = .5) +
  geom_hline(yintercept = 1.936591) +
  ylab("Child RSV B VARscore") +
  theme_light()+
  scale_color_manual(values = c("black", "#005A32", "#41AB5D", "#A1D99B")) +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())

wilcox.test((hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
               left_join(Maternal_pos) %>% filter(taxon_species == "Human respiratory syncytial virus B") %>%
               filter(!is.na(`Mom RSVA quartile`)) %>% filter(`Mom RSVA quartile` == "1st quartile"))$vir_score,
            (hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
               left_join(Maternal_pos) %>% filter(taxon_species == "Human respiratory syncytial virus B") %>%
               filter(!is.na(`Mom RSVA quartile`)) %>% filter(`Mom RSVA quartile` == "4th quartile"))$vir_score) # not significantly different

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - RSV B VARscore dotplot by maternal RSV A antibodies.jpeg",
       dotplot,
       scale = 1,
       height = 2,
       width = 3.5,
       units = "in",
       dpi = 300)


positive_rate <- hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
  left_join(Maternal_pos) %>% filter(taxon_species == "Human respiratory syncytial virus B") %>%
  filter(!is.na(`Mom RSVA quartile`)) %>%
  mutate(positive = vir_score > 1.936591) %>%
  group_by(`Mom RSVA quartile`) %>%
  summarise(n_positive = sum(positive), n_total = n()) %>%
  mutate(proportion = n_positive / n_total)

fisher.test(matrix(data = c(110, 132,
                            109, 132), nrow = 2))
