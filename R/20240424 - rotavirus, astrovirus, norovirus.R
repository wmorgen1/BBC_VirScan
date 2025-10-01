setwd("~/R/Larman_Lab")
library(tidyverse)
library(RColorBrewer)
library(gtable)
library(grid)
library(ggplot2)
library(ggbeeswarm)
source("ArboScan/common_functions.R")
library(data.table)
#library(ARscore)

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
  filter(taxon_species %in% c("Rotavirus A",
                              "Rotavirus X",
                              "Human astrovirus-1",
                              "Human astrovirus-8",
                              "Norwalk virus",
                              "Lordsdale virus")) %>%
  mutate(taxon_species = as.factor(taxon_species)) %>%
  arrange(taxon_species)

dotplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>% group_by(Plasma_type),
             aes(x = interaction(Plasma_type, taxon_species), y = vir_score,
                 color = Plasma_type),
             position = position_quasirandom(),
             size = .1,
             alpha = .3) +
  ylab("VARscore") +
  xlab("Virus") +
  geom_hline(yintercept = 1.936591) +
  scale_x_discrete(labels=c("cord.Rotavirus A" = "Rota A",
                            "child.Rotavirus A" = "Rota A",
                            "cord.Rotavirus X" = "Rota X",
                            "child.Rotavirus X" = "Rota X",
                            "cord.Human astrovirus-1" = "Astro 1",
                            "child.Human astrovirus-1" = "Astro 1",
                            "cord.Human astrovirus-8" = "Astro 8",
                            "child.Human astrovirus-8" = "Astro 8",
                            "cord.Norwalk virus" = "Norwalk",
                            "child.Norwalk virus" = "Norwalk",
                            "cord.Lordsdale virus" = "Lordsdale",
                            "child.Lordsdale virus" = "Lordsdale"
                            )) +
  theme_light() +
  scale_color_manual(values = c("black", "red")) +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - rota noro astro VARscore dotplot.jpeg",
       dotplot,
       scale = 1,
       height = 2,
       width = 7.5,
       units = "in",
       dpi = 300)

##

timeplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% filter(Plasma_type == "child") %>%
               filter(taxon_species %in% c("Rotavirus A",
                                           "Human astrovirus-1",
                                           "Norwalk virus", 
                                           "Lordsdale virus")) %>%
               mutate(taxon_species = case_when(
                 taxon_species == "Rotavirus A" ~ "Rota A", 
                 taxon_species == "Human astrovirus-1" ~ "Astro 1",
                 taxon_species == "Norwalk virus" ~ "Norwalk",
                 taxon_species == "Lordsdale virus" ~ "Lordsdale"
               )),
             aes(x = years, y = vir_score),
             size = .3,
             alpha = .3) +
  geom_hline(yintercept = 1.936591) +
  ylab("VARscore") +
  xlab("Age (years)") +
  theme_light() + 
  facet_wrap(vars(taxon_species), ncol = 2) + 
  theme(text = element_text(size = 8))

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - child rota astro noro virus VARscore vs time dotplot.jpeg",
       timeplot,
       scale = 1,
       height = 3,
       width = 3,
       units = "in",
       dpi = 300)

hsv_serobyage = hsv_VAR %>% filter(Plasma_type == "child") %>%
  mutate(age_group = case_when(
    years < 1 ~ "0 - 1",
    years < 2 ~ "1 - 2",
    years < 3 ~ "2 - 3",
    years < 4 ~ "3 - 4",
    years >= 4 ~ "4+"
  )) %>%
  mutate(seropositive = (vir_score > 1.936591)) %>%
  ungroup() %>% group_by(age_group, taxon_species) %>%
  summarise(proportion = sum(seropositive)/n(), number = n())

sero_byage <- ggplot() + 
  geom_col(data = hsv_serobyage,
             aes(x = taxon_species, y = proportion,
                 fill = age_group),
             position = position_dodge(),
             color = "black"
             ) +
  ylab("Proportion positive") +
  xlab("Age (years)") +
  theme_light() + 
  theme(text = element_text(size = 12)) +
  scale_x_discrete(labels=c("Rotavirus A" = "Rota A",
                            "Rotavirus X" = "Rota X",
                            "Human astrovirus-1" = "Astro 1",
                            "Human astrovirus-8" = "Astro 8",
                            "Norwalk virus" = "Norwalk",
                            "Lordsdale virus" = "Lordsdale"
  )) +
  scale_fill_viridis_d() +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - child rota astro noro seropositivity by age.jpeg",
       sero_byage,
       scale = 1,
       height = 1.5,
       width = 3.5,
       units = "in",
       dpi = 300)

##############
# cross reactivity in kids
##############

HSV1_HSV2 <- hsv_VAR %>% ungroup() %>%
  pivot_wider(id_cols = c(ID, Plasma_type), names_from = taxon_species,
              values_from = vir_score, values_fn = max)

HSV_xreact <- ggplot() + 
  geom_point(data = HSV1_HSV2,
             aes(x = `Lordsdale virus`, 
                 y = `Norwalk virus`,
                 color = Plasma_type),
             size = .3,
             alpha = .3) +
  xlab("Lordsdale") +
  ylab("Norwalk") +
  theme_light() +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(vars(Plasma_type)) +
  theme(text = element_text(size = 8),
        legend.position = "none") 

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - Noro cross reactivity.jpeg",
       HSV_xreact,
       scale = 1,
       height = 1.5,
       width = 3,
       units = "in",
       dpi = 300)

library(GGally)


ggpairs(data = HSV1_HSV2,
        columns = c(3:8),
        #mapping = aes(alpha = .1)
        lower = list(continuous = pairs_dots)) +
  theme_light()
