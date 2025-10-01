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

# Enterovirus
entero_max_VAR <- BBC_VARscores %>%
  filter(taxon_genus %in% c("Enterovirus")) %>%
  group_by(taxon_genus, ID, Plasma_type, lab_id, years) %>%
  summarise(vir_score = max(vir_score))

# parechovirus
parecho_max_VAR <- BBC_VARscores %>%
  filter(taxon_genus %in% c("Parechovirus")) %>%
  group_by(taxon_genus, ID, Plasma_type, lab_id, years) %>%
  summarise(vir_score = max(vir_score))

# Mastadenovirus
adeno_max_VAR <- BBC_VARscores %>%
  filter(taxon_species %in% c("Human adenovirus A serotype 12", 
                              "Human adenovirus B serotype 16",
                              "Human adenovirus C serotype 2",
                              "Human adenovirus D37",
                              "Human adenovirus E serotype 4", 
                              "Human adenovirus F serotype 40")) %>%
  group_by(taxon_genus, ID, Plasma_type, lab_id, years) %>%
  summarise(vir_score = max(vir_score))

# Norovirus
noro_max_VAR <- BBC_VARscores %>%
  filter(taxon_genus %in% c("Norovirus")) %>%
  group_by(taxon_genus, ID, Plasma_type, lab_id, years) %>%
  summarise(vir_score = max(vir_score))

# Astrovirus
astro_max_VAR <- BBC_VARscores %>%
  filter(taxon_genus %in% c("Mamastrovirus")) %>%
  group_by(taxon_genus, ID, Plasma_type, lab_id, years) %>%
  summarise(vir_score = max(vir_score))

# Rotavirus
rota_max_VAR <- BBC_VARscores %>%
  filter(taxon_genus %in% c("Rotavirus")) %>%
  group_by(taxon_genus, ID, Plasma_type, lab_id, years) %>%
  summarise(vir_score = max(vir_score))

hsv_VAR <- entero_max_VAR %>%
  full_join(parecho_max_VAR) %>%
  full_join(adeno_max_VAR) %>%
  full_join(noro_max_VAR) %>%
  full_join(astro_max_VAR) %>%
  full_join(rota_max_VAR) %>%
  mutate(taxon_genus = as.factor(taxon_genus)) %>%
  mutate(taxon_genus = fct_relevel(taxon_genus,
                                   c("Enterovirus",
                                   "Parechovirus",
                                   "Mastadenovirus",
                                   "Norovirus" ,
                                   "Mamastrovirus",
                                   "Rotavirus")))

dotplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>% group_by(Plasma_type),
             aes(x = interaction(Plasma_type, taxon_genus), y = vir_score,
                 color = Plasma_type),
             position = position_quasirandom(),
             size = .1,
             alpha = .3) +
  ylab("Max VARscore in Genus") +
  xlab("Virus") +
  geom_hline(yintercept = 1.936591) +
  scale_x_discrete(labels=c("cord.Enterovirus" = "EV",
                            "child.Enterovirus" = "EV",
                            "cord.Parechovirus" = "pEV",
                            "child.Parechovirus" = "pEV",
                            "cord.Mastadenovirus" = "Adeno",
                            "child.Mastadenovirus" = "Adeno",
                            "cord.Norovirus" = "Noro",
                            "child.Norovirus" = "Noro",
                            "cord.Mamastrovirus" = "Astro",
                            "child.Mamastrovirus" = "Astro",
                            "cord.Rotavirus" = "Rota",
                            "child.Rotavirus" = "Rota"
                            ),
                   breaks = c("child.Enterovirus",
                              "cord.Enterovirus",
                              "child.Parechovirus",
                              "cord.Parechovirus",
                              "child.Mastadenovirus",
                              "cord.Mastadenovirus" ,
                              "child.Norovirus",
                              "cord.Norovirus",
                              "child.Mamastrovirus",
                              "cord.Mamastrovirus",
                              "child.Rotavirus",
                              "cord.Rotavirus")) +
  theme_light() +
scale_color_manual(values = c("black", "red")) +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.position = "none")

ggsave("Figures/20240515 - Fig 3 genus VARscore dotplot.jpeg",
       dotplot,
       scale = 1,
       height = 2,
       width = 6.5,
       units = "in",
       dpi = 300)

seropositivity <- hsv_VAR %>% ungroup() %>% group_by(taxon_genus, Plasma_type) %>%
  mutate(positive = (vir_score > 1.936591)) %>%
  summarise(n = n(), positive = sum(positive)) %>%
  mutate(percent = positive / n)
  

timeplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>%
               filter(Plasma_type == "child") %>%
               mutate(taxon_genus = case_when(
                 taxon_genus == "Enterovirus" ~ "EV",
                 taxon_genus == "Parechovirus" ~ "pEV",
                 taxon_genus == "Mastadenovirus" ~ "Adeno",
                 taxon_genus == "Norovirus" ~ "Noro",
                 taxon_genus == "Mamastrovirus" ~ "Astro",
                 taxon_genus == "Rotavirus" ~ "Rota"
               )) %>% 
               mutate(taxon_genus = as.factor(taxon_genus)) %>%
               mutate(taxon_genus = fct_relevel(taxon_genus,
                                                c("EV",
                                                  "pEV",
                                                  "Adeno",
                                                  "Noro" ,
                                                  "Astro",
                                                  "Rota"))),
             aes(x = years, y = vir_score),
             size = .3,
             alpha = .3) +
  geom_hline(yintercept = 1.936591) +
  ylab("VARscore") +
  xlab("Age (years)") +
  theme_light() + 
  facet_wrap(vars(taxon_genus), nrow = 1) + 
  theme(text = element_text(size = 8))

ggsave("Figures/20240515 - child genus VARscore vs time dotplot.jpeg",
       timeplot,
       scale = 1,
       height = 1.5,
       width = 6.5,
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

sero_byage <- ggplot() + 
  geom_col(data = hsv_serobyage %>%
             mutate(taxon_genus = case_when(
               taxon_genus == "Enterovirus" ~ "EV",
               taxon_genus == "Parechovirus" ~ "pEV",
               taxon_genus == "Mastadenovirus" ~ "Adeno",
               taxon_genus == "Norovirus" ~ "Noro",
               taxon_genus == "Mamastrovirus" ~ "Astro",
               taxon_genus == "Rotavirus" ~ "Rota"
             )) %>% 
             mutate(taxon_genus = as.factor(taxon_genus)) %>%
             mutate(taxon_genus = fct_relevel(taxon_genus,
                                              c("EV",
                                                "pEV",
                                                "Adeno",
                                                "Noro" ,
                                                "Astro",
                                                "Rota"))),
           aes(x = taxon_genus, y = proportion,
               fill = age_group),
           position = position_dodge(),
             color = "black"
             ) +
  #scale_x_discrete(labels=c("Enterovirus A" = "EVA",
  #                          "Enterovirus B" = "EVB",
  #                          "Enterovirus C" = "EVC",
  #                          "Enterovirus D" = "EVD",
  #                          "Rhinovirus A" = "RVA",
  #                          "Rhinovirus B" = "RVB",
  #                          "Human parechovirus 1" = "pEV1",
  #                          "Human parechovirus 2" = "pEV2"
  #)) +
  ylab("Proportion positive") +
  theme_light() + 
  scale_fill_viridis_d() +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240515 - child genus seropositivity by age.jpeg",
       sero_byage,
       scale = 1,
       height = 1.5,
       width = 4,
       units = "in",
       dpi = 300)

##############
# cross reactivity in kids
##############

HSV1_HSV2 <- hsv_VAR %>% ungroup() %>%
  filter(taxon_genus == "Enterovirus") %>%
  pivot_wider(id_cols = c(ID, Plasma_type), names_from = taxon_species,
              values_from = vir_score, values_fn = max)

xreact <- ggplot() + 
  geom_point(data = HSV1_HSV2,
             aes(x = `Enterovirus A`, 
                 y = `Enterovirus B`,
                 color = Plasma_type),
             size = .3,
             alpha = .3) +
  xlab("EVA") +
  ylab("EVB") +
  theme_light() +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(vars(Plasma_type), ncol = 1) +
  theme(text = element_text(size = 8),
        legend.position = "none")

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - EVA EVB cross reactivity.jpeg",
       xreact,
       scale = 1,
       height = 3,
       width = 2,
       units = "in",
       dpi = 300)

xreact_2 <- ggplot() + 
  geom_point(data = HSV1_HSV2,
             aes(x = `Enterovirus A`, 
                 y = `Enterovirus B`,
                 color = Plasma_type),
             size = .3,
             alpha = .3) +
  xlab("EVA") +
  ylab("EVB") +
  theme_light() +
  scale_color_manual(values = c("red", "black")) +
  facet_wrap(vars(Plasma_type), ncol = 1) +
  theme(text = element_text(size = 8),
        legend.position = "none")
