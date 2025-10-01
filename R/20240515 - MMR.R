setwd("C:/Users/wmorgen1/Documents/Larman_Lab/BBC/")
library(tidyverse)
library(RColorBrewer)
library(gtable)
library(grid)
library(ggplot2)
library(ggbeeswarm)
#source("ArboScan/common_functions.R")
library(data.table)
#library(ARscore)

########
# Data Prep
#######
library(readxl)
sample_meta <- read_excel("sample_meta.xls")

sample_meta <- sample_meta %>% mutate(age = difftime(time1 = Plasma_Date, time2 = DOD, units = "days")) %>%
  mutate(years = (age %>% as.numeric())/365.25) %>% mutate(Plasma_type = ifelse(Plasma_type=="Postnatal","child",Plasma_type))

BBC_VARscores <- read_rds("Data/20240424 - BBC VARscores - maternal abs removed.rds")

BBC_VARscores <- BBC_VARscores %>% dplyr::rename(vir_score = value) %>%
  left_join(sample_meta, multiple = "all")

hsv_VAR <- BBC_VARscores %>%
  filter(taxon_genus %in% c("Morbillivirus", "Rubulavirus", "Rubivirus")) %>%
  filter(taxon_species %in% c("Measles virus", "Rubella virus",
                              "Mumps virus")) %>%
  mutate(taxon_species = as.factor(taxon_species)) %>%
  mutate(taxon_species = fct_relevel(taxon_species,
                                     "Measles virus", "Rubella virus",
                                     "Mumps virus")) %>%
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
  scale_x_discrete(labels=c("cord.Measles virus" = "Measles",
                            "child.Measles virus" = "Measles",
                            "cord.Rubella virus" = "Rubella",
                            "child.Rubella virus" = "Rubella",
                            "cord.Mumps virus" = "Mumps",
                            "child.Mumps virus" = "Mumps")) +
  theme_light() +
  scale_color_manual(values = c("black", "red")) +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240515 - MMR VARscore dotplot.jpeg",
       dotplot,
       scale = 1,
       height = 2,
       width = 4,
       units = "in",
       dpi = 300)

seropositivity <- hsv_VAR %>% ungroup() %>% group_by(taxon_genus, Plasma_type) %>%
  mutate(positive = (vir_score > 1.936591)) %>%
  summarise(n = n(), positive = sum(positive)) %>%
  mutate(percent = positive / n)


timeplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% filter(Plasma_type == "child") %>%
               mutate(taxon_species = case_when(
                 taxon_species == "Measles virus" ~ "Measles",
                 taxon_species == "Rubella virus" ~ "Rubella",
                 taxon_species == "Mumps virus" ~ "Mumps")),
             aes(x = years, y = vir_score),
             size = .3,
             alpha = .3,) +
  geom_hline(yintercept = 1.936591) +
  ylab("VARscore") +
  xlab("Age (years)") +
  theme_light() + 
  facet_wrap(vars(taxon_species)) + 
  theme(text = element_text(size = 8)) 

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240515 - child MMR VARscore vs time dotplot.jpeg",
       timeplot,
       scale = 1,
       height = 1.5,
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
  ungroup() %>% group_by(age_group, taxon_genus) %>%
  summarise(proportion = sum(seropositive)/n(), positive = sum(seropositive), number = n())

mmr_serobyage = hsv_VAR %>% filter(Plasma_type == "child") %>%
  filter(taxon_species != "Human parainfluenza 4a virus") %>%
  group_by(ID, age, years) %>% summarise(vir_score = max(vir_score)) %>%
  ungroup() %>% 
  mutate(age_group = case_when(
    years < 1 ~ "0 - 1",
    years < 2 ~ "1 - 2",
    years < 3 ~ "2 - 3",
    years < 4 ~ "3 - 4",
    years >= 4 ~ "4+"
  )) %>%
  mutate(seropositive = (vir_score > 1.936591)) %>%
  ungroup() %>% group_by(age_group) %>%
  summarise(proportion = sum(seropositive)/n(), positive = sum(seropositive), number = n())

mmr_sero = hsv_VAR %>% 
  filter(taxon_species != "Human parainfluenza 4a virus") %>%
  group_by(ID, age, years, Plasma_type) %>% summarise(vir_score = max(vir_score)) %>%
  ungroup() %>%
  mutate(seropositive = (vir_score > 1.936591)) %>%
  ungroup() %>% group_by(Plasma_type) %>%
  summarise(proportion = sum(seropositive)/n(), positive = sum(seropositive), number = n())

sero_byage <- ggplot() + 
  geom_col(data = mmr_serobyage,
           aes(x = '1', y = proportion,
               fill = age_group),
           position = position_dodge(),
           color = "black"
  ) +
  ylab("Proportion positive") +
  scale_x_discrete(labels=c("1" = "Any MMR"
  )) +
  theme_light() + 
  scale_fill_viridis_d() +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240515 - MMR seropositivity by age.jpeg",
       sero_byage,
       scale = 1,
       height = 1.5,
       width = 2,
       units = "in",
       dpi = 300)

##############
# cross reactivity in kids
##############

HSV1_HSV2 <- hsv_VAR %>% ungroup() %>%
  pivot_wider(id_cols = c(ID, Plasma_type), names_from = taxon_species,
              values_from = vir_score, values_fn = max)

HSV_xreact <- ggplot() + 
  geom_point(data = HSV1_HSV2 %>% filter(Plasma_type == "child"),
             aes(x = `Measles virus`, 
                 y = `Rubella virus`),
             size = .3,
             alpha = .3) +
  xlab("Measles") +
  ylab("Rubella") +
  theme_light() +
  facet_wrap(vars(Plasma_type)) +
  theme(text = element_text(size = 8),
        legend.position = "none") #+ 
 #geom_hline(yintercept = 1) + 
 #geom_vline(xintercept = 1)

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240515 - MMR co-exposure.jpeg",
       HSV_xreact,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300)

cor.test((HSV1_HSV2 %>% filter(Plasma_type == "child"))$`Measles virus`,
         (HSV1_HSV2 %>% filter(Plasma_type == "child"))$`Rubella virus`)
# pearson's r = 0.68, p < 2.2e-16
