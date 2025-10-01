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
  filter(taxon_genus %in% c("Influenzavirus A",
                            "Influenzavirus B",
                            "Influenzavirus C")) %>%
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
  scale_x_discrete(labels=c("cord.Influenza A virus" = "Flu A",
                            "child.Influenza A virus" = "Flu A",
                            "cord.Influenza B virus" = "Flu B",
                            "child.Influenza B virus" = "Flu B",
                            "cord.Influenza C virus" = "Flu C",
                            "child.Influenza C virus" = "Flu C"
                            )) +
  theme_light() +
  scale_color_manual(values = c("black", "red")) +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())
#,       axis.text.x = element_text(angle = 90))

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - influenza VARscore dotplot.jpeg",
       dotplot,
       scale = 1,
       height = 2,
       width = 3.5,
       units = "in",
       dpi = 300)

seropositivity <- hsv_VAR %>% ungroup() %>% group_by(taxon_genus, Plasma_type) %>%
  mutate(positive = (vir_score > 1.936591)) %>%
  summarise(n = n(), positive = sum(positive)) %>%
  mutate(percent = positive / n)


timeplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% filter(Plasma_type == "child") %>%
               mutate(taxon_species = case_when(
                 taxon_species == "Influenza A virus" ~ "Flu A",
                 taxon_species == "Influenza B virus" ~ "Flu B",
                 taxon_species == "Influenza C virus" ~ "Flu C"
               )),
             aes(x = years, y = vir_score),
             size = .3,
             alpha = .3) +
  geom_hline(yintercept = 1.936591) +
  ylab("VARscore") +
  xlab("Age (years)") +
  theme_light() + 
  facet_wrap(vars(taxon_species)) + 
  theme(text = element_text(size = 8))

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - child Flu VARscore vs time dotplot.jpeg",
       timeplot,
       scale = 1,
       height = 1.5,
       width = 4,
       units = "in",
       dpi = 300)

dateplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% filter(Plasma_type == "child") %>%
               mutate(taxon_species = case_when(
                 taxon_species == "Influenza A virus" ~ "Flu A",
                 taxon_species == "Influenza B virus" ~ "Flu B",
                 taxon_species == "Influenza C virus" ~ "Flu C"
               )),
             aes(x = `Plasma_Date`, y = vir_score),
             size = 1,
             alpha = .3) +
  ylab("VARscore") +
  xlab("Date") +
  theme_light() + 
  geom_smooth(data = hsv_VAR %>% filter(Plasma_type == "child") %>%
                mutate(taxon_species = case_when(
                  taxon_species == "Influenza A virus" ~ "Flu A",
                  taxon_species == "Influenza B virus" ~ "Flu B",
                  taxon_species == "Influenza C virus" ~ "Flu C"
                )),
              aes(x = `Plasma_Date`, y = vir_score),
              formula = y ~ x,
              span = 0.25) +
  facet_wrap(vars(taxon_species))

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - child Flu VARscore vs date dotplot - seasonality and outbreaks.jpeg",
       dateplot,
       scale = 1,
       height = 2,
       width = 4,
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
  summarise(proportion = sum(seropositive)/n(), positive = sum(seropositive), number = n())

sero_byage <-ggplot() + 
  geom_col(data = hsv_serobyage,
             aes(x = taxon_species, y = proportion,
                 fill = age_group),
             position = position_dodge(),
             color = "black"
             ) +
  ylab("Proportion positive") +
  xlab("Age (years)") +
  scale_x_discrete(labels=c("Influenza A virus" = "Flu A",
                            "Influenza B virus" = "Flu B",
                            "Influenza C virus" = "Flu C")) +
  theme_light() + 
  scale_fill_viridis_d() +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.position = "none")

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - child Influenza seropositivity by age.jpeg",
       sero_byage,
       scale = 1,
       height = 1,
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
  geom_point(data = HSV1_HSV2,
             aes(x = `Influenza A virus`, 
                 y = `Influenza B virus`,
                 color = Plasma_type),
             size = .3,
             alpha = .3) +
  xlab("Flu A") +
  ylab("Flu B") +
  theme_light() +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(vars(Plasma_type)) +
  theme(text = element_text(size = 8),
        legend.position = "none")

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - Flu A B cross reactivity.jpeg",
       HSV_xreact,
       scale = 1,
       height = 1.5,
       width = 3,
       units = "in",
       dpi = 300)
