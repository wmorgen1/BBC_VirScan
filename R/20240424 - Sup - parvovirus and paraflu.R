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
   filter(taxon_genus == "Erythroparvovirus") %>%
  mutate(taxon_species = as.factor(taxon_species)) %>%
  arrange(taxon_species)

ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>% group_by(Plasma_type),
             aes(x = interaction(Plasma_type, taxon_species), y = vir_score,
                 color = Plasma_type),
             position = position_quasirandom(),
             size = 1,
             alpha = .3) +
  ylab("VARscore") +
  xlab("Virus") +
  scale_x_discrete(labels=c("cord.Human parvovirus B19" = "B19",
                            "child.Human parvovirus B19" = "B19"
                            )) +
  theme_light()

seropositivity <- hsv_VAR %>% ungroup() %>% group_by(taxon_genus, Plasma_type) %>%
  mutate(positive = (vir_score > 1.936591)) %>%
  summarise(n = n(), positive = sum(positive)) %>%
  mutate(percent = positive / n)

ggplot() + 
  geom_point(data = hsv_VAR %>% filter(Plasma_type == "child"),
             aes(x = years, y = vir_score),
             size = 1,
             alpha = .3) +
  ylab("VARscore") +
  xlab("Age (years)") +
  theme_light() + 
  #geom_smooth(data = hsv_VAR %>% filter(Plasma_type == "Postnatal"),
  #            aes(x = years, y = vir_score),
  #            method = "lm",
  #            formula = y ~ x) +
  facet_wrap(vars(taxon_species)) + 
  theme(text = element_text(size = 12))

ggplot() + 
  geom_point(data = hsv_VAR %>% filter(Plasma_type == "Postnatal"),
             aes(x = `Plasma_Date`, y = vir_score),
             size = 1,
             alpha = .3) +
  ylab("VARscore") +
  xlab("Date") +
  theme_light() + 
  #geom_smooth(data = hsv_VAR %>% filter(Plasma_type == "Postnatal"),
  #            aes(x = `Plasma_Date`, y = vir_score),
  #            method = "lm",
  #            formula = y ~ x) +
  facet_wrap(vars(taxon_species))

hsv_serobyage = hsv_VAR %>% filter(Plasma_type == "child") %>%
  mutate(age_group = case_when(
    years < 1 ~ "0 - 1",
    years < 2 ~ "1 - 2",
    years < 3 ~ "2 - 3",
    years < 4 ~ "3 - 4",
    years >= 4 ~ "4+"
  )) %>%
  mutate(seropositive = (vir_score > 2)) %>%
  ungroup() %>% group_by(age_group, taxon_species) %>%
  summarise(proportion = sum(seropositive)/n(),positive = sum(seropositive), number = n())

######
hsv_VAR <- BBC_VARscores %>%
  filter(taxon_species %in% c("Human parainfluenza 1 virus",
                              "Human parainfluenza 2 virus",
                              "Human parainfluenza 3 virus",
                              "Human parainfluenza 4a virus")) %>%
  mutate(taxon_genus = "pf") %>%
  group_by(taxon_genus, ID, Plasma_type, lab_id, years) %>%
  summarise(vir_score = max(vir_score))


dotplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>% group_by(Plasma_type),
             aes(x = interaction(Plasma_type, taxon_genus), y = vir_score,
                 color = Plasma_type),
             position = position_quasirandom(),
             size = .1,
             alpha = .3) +
  ylab("Max VARscore") +
  xlab("Virus") +
  geom_hline(yintercept = 1.936591) +
  scale_x_discrete(labels=c("cord.pf" = "pFlu",
                            "child.pf" = "pFlu"
  ),
  breaks = c("child.pf",
             "cord.pf")) +
  theme_light() +
  scale_color_manual(values = c("black", "red")) +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.position = "none")

ggsave("Figures/20240515 - Fig S7 paraflu VARscore dotplot.jpeg",
       dotplot,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300)

seropositivity <- hsv_VAR %>% ungroup() %>% group_by(taxon_genus, Plasma_type) %>%
  mutate(positive = (vir_score > 1.936591)) %>%
  summarise(n = n(), positive = sum(positive)) %>%
  mutate(percent = positive / n)

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
