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

# excluding unspecified
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
  geom_point(data = hsv_VAR %>% ungroup() %>% group_by(Plasma_type),
             aes(x = interaction(Plasma_type, taxon_species), y = vir_score,
                 color = Plasma_type),
             position = position_quasirandom(),
             size = .1,
             alpha = .3) +
  geom_hline(yintercept = 1.936591) +
  ylab("VARscore") +
  xlab("Virus") +
  scale_x_discrete(labels=c("cord.Human metapneumovirus" = "MPnV",
                            "child.Human metapneumovirus" = "MPnV",
                            "cord.Human respiratory syncytial virus A" = "RSV A",
                            "child.Human respiratory syncytial virus A" = "RSV A",
                            "cord.Human respiratory syncytial virus B" = "RSV B",
                            "child.Human respiratory syncytial virus B" = "RSV B"
                            )) +
  theme_light()+
  scale_color_manual(values = c("black", "red")) +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())
#,       axis.text.x = element_text(angle = 90))

seropositivity <- hsv_VAR %>% ungroup() %>%
  group_by(taxon_genus, Plasma_type, ID) %>%
  summarise(vir_score = max(vir_score)) %>%
  mutate(positive = (vir_score > 1.936591)) %>%
  summarise(n = n(), positive = sum(positive)) %>%
  mutate(percent = positive / n)

seropositivity <- hsv_VAR %>% ungroup() %>%
  group_by(taxon_species, Plasma_type) %>%
  mutate(positive = (vir_score > 1.936591)) %>%
  summarise(n = n(), positive = sum(positive)) %>%
  mutate(percent = positive / n)

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - pneumovirus VARscore dotplot.jpeg",
       dotplot,
       scale = 1,
       height = 2,
       width = 3.5,
       units = "in",
       dpi = 300)

timeplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% filter(Plasma_type == "child") %>%
               mutate(taxon_species = case_when(
                 taxon_species == "Human metapneumovirus" ~ "MPnV",
                 taxon_species == "Human respiratory syncytial virus A" ~ "RSV A",
                 taxon_species == "Human respiratory syncytial virus B" ~ "RSV B"
               )),
             aes(x = years, y = vir_score),
             size = .3,
             alpha = .3,
             color = "black") +
  geom_hline(yintercept = 1.936591) +
  ylab("VARscore") +
  xlab("Age (years)") +
  theme_light() + 
  #geom_smooth(data = hsv_VAR %>% filter(Plasma_type == "child"),
  #            aes(x = years, y = vir_score),
  #            method = "lm",
  #            formula = y ~ x) +
  facet_wrap(vars(taxon_species)) + 
  theme(text = element_text(size = 8))

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - child pneumovirus VARscore vs time dotplot.jpeg",
       timeplot,
       scale = 1,
       height = 1.5,
       width = 4,
       units = "in",
       dpi = 300)

dateplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% filter(Plasma_type == "child")%>%
               mutate(taxon_species = case_when(
                 taxon_species == "Human metapneumovirus" ~ "MPnV",
                 taxon_species == "Human respiratory syncytial virus A" ~ "RSV A",
                 taxon_species == "Human respiratory syncytial virus B" ~ "RSV B"
               )),
             aes(x = `Plasma_Date`, y = vir_score),
             size = 1,
             alpha = .3) +
  ylab("VARscore") +
  xlab("Date") +
  theme_light() + 
  geom_smooth(data = hsv_VAR %>% filter(Plasma_type == "child")%>%
                mutate(taxon_species = case_when(
                  taxon_species == "Human metapneumovirus" ~ "MPnV",
                  taxon_species == "Human respiratory syncytial virus A" ~ "RSV A",
                  taxon_species == "Human respiratory syncytial virus B" ~ "RSV B"
                )),
              aes(x = `Plasma_Date`, y = vir_score),
              formula = y ~ x,
              span = 0.25) +
  facet_wrap(vars(taxon_species))

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - child pneumo VARscore vs date dotplot - seasonality and outbreaks.jpeg",
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

sero_byage <- ggplot() + 
  geom_col(data = hsv_serobyage,
             aes(x = taxon_species, y = proportion,
                 fill = age_group),
             position = position_dodge(),
             color = "black"
             ) +
  ylab("Proportion positive") +
  xlab("Age (years)") +
  scale_x_discrete(labels=c("Human metapneumovirus" = "MPnV",
                            "Human respiratory syncytial virus A" = "RSV A",
                            "Human respiratory syncytial virus B" = "RSV B"
  )) +
  theme_light() + 
  scale_fill_viridis_d() +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.position = "none")

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - child pneumovirus seropositivity by age.jpeg",
       sero_byage,
       scale = 1,
       height = 1,
       width = 2,
       units = "in",
       dpi = 300)

hsv_serobyage = hsv_VAR %>% filter(Plasma_type == "child") %>%
  group_by(taxon_genus, Plasma_type, years, ID) %>%
  summarise(vir_score = max(vir_score)) %>%
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

##############
# cross reactivity 
##############

HSV1_HSV2 <- hsv_VAR %>% ungroup() %>%
  pivot_wider(id_cols = c(ID, Plasma_type), names_from = taxon_species,
              values_from = vir_score, values_fn = max)

HSV_xreact <- ggplot() + 
  geom_point(data = HSV1_HSV2,
             aes(x = `Human respiratory syncytial virus A`, 
                 y = `Human respiratory syncytial virus B`,
                 color = Plasma_type),
             size = .3,
             alpha = .3) +
  xlab("RSV A") +
  ylab("RSV B") +
  theme_light() +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(vars(Plasma_type)) +
  theme(text = element_text(size = 8),
        legend.position = "none") #+ 
  #geom_hline(yintercept = 1.936591)+ 
  #geom_vline(xintercept = 1.936591)+
  #geom_abline(intercept = 0, slope = 2)+
  #geom_abline(intercept = 0, slope = 0.5)


ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - RSV cross reactivity.jpeg",
       HSV_xreact,
       scale = 1,
       height = 1.5,
       width = 3,
       units = "in",
       dpi = 300)
