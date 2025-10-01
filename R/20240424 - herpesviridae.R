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
  filter(taxon_species %in% c("Human herpesvirus 1", "Human herpesvirus 2",
                              "Human herpesvirus 3",
                              "Human cytomegalovirus", "Epstein-Barr virus",
                              "Human herpesvirus 6A", "Human herpesvirus 6B",
                              "Human herpesvirus 7", "Human herpesvirus 8")) %>%
  mutate(taxon_species = as.factor(taxon_species)) %>%
  mutate(taxon_species = fct_relevel(taxon_species, "Human herpesvirus 1", "Human herpesvirus 2",
                                     "Human herpesvirus 3",
                                     "Epstein-Barr virus", "Human cytomegalovirus",
                                     "Human herpesvirus 6A", "Human herpesvirus 6B",
                                     "Human herpesvirus 7", "Human herpesvirus 8")) %>%
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
  scale_x_discrete(labels=c("cord.Human herpesvirus 1" = "HSV1",
                            "child.Human herpesvirus 1" = "HSV1",
                            "cord.Human herpesvirus 2" = "HSV2",
                            "child.Human herpesvirus 2" = "HSV2",
                            #"cord.Varicella-zoster virus" = "VZV",
                            #"child.Varicella-zoster virus" = "VZV",
                            "cord.Human herpesvirus 3" = "HHV3",
                            "child.Human herpesvirus 3" = "HHV3",
                            "cord.Epstein-Barr virus" = "EBV",
                            "child.Epstein-Barr virus" = "EBV",
                            "cord.Human cytomegalovirus" = "CMV",
                            "child.Human cytomegalovirus" = "CMV",
                            "cord.Human herpesvirus 6A" = "HHV6A",
                            "child.Human herpesvirus 6A" = "HHV6A",
                            "cord.Human herpesvirus 6B" = "HHV6B",
                            "child.Human herpesvirus 6B" = "HHV6B",
                            "cord.Human herpesvirus 7" = "HHV7",
                            "child.Human herpesvirus 7" = "HHV7",
                            "cord.Human herpesvirus 8" = "HHV8",
                            "child.Human herpesvirus 8" = "HHV8"
                            )) +
  theme_light() +
  scale_color_manual(values = c("black", "red")) +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.position = "none")

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - herpesvirus VARscore dotplot.jpeg",
       dotplot,
       scale = 1,
       height = 2,
       width = 7,
       units = "in",
       dpi = 300)

seropositivity <- hsv_VAR %>% ungroup() %>% group_by(taxon_species, Plasma_type) %>%
  mutate(positive = (vir_score > 1.936591)) %>%
  summarise(n = n(), positive = sum(positive)) %>%
  mutate(percent = positive / n)


timeplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% filter(Plasma_type == "child") %>%
               filter(taxon_species %in% c("Human herpesvirus 1", "Human herpesvirus 2",
                                           "Human herpesvirus 3",
                                           "Human cytomegalovirus", "Epstein-Barr virus",
                                           "Human herpesvirus 6B")) %>%
               mutate(taxon_species = case_when(
                 taxon_species == "Human herpesvirus 1" ~ "HSV 1", 
                 taxon_species == "Human herpesvirus 2" ~ "HSV 2",
                 taxon_species == "Human herpesvirus 3" ~ "VZV",
                 taxon_species == "Human cytomegalovirus" ~ "CMV", 
                 taxon_species == "Epstein-Barr virus" ~ "EBV",
                 taxon_species == "Human herpesvirus 6B" ~ "HHV6B"
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

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - child herpesvirus VARscore vs time dotplot.jpeg",
       timeplot,
       scale = 1,
       height = 3,
       width = 4,
       units = "in",
       dpi = 300)

# check for seasonality, outbreaks, timing with vaccination...
ggplot() + 
  geom_point(data = hsv_VAR %>% filter(Plasma_type == "child"),
             aes(x = `Plasma_Date`, y = vir_score),
             size = 1,
             alpha = .3) +
  ylab("VARscore") +
  xlab("Date") +
  theme_light() + 
  #geom_smooth(data = hsv_VAR %>% filter(Plasma_type == "child"),
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
  mutate(seropositive = (vir_score > 1.936591)) %>%
  ungroup() %>% group_by(age_group, taxon_species) %>%
  summarise(proportion = sum(seropositive)/n(), positive = sum(seropositive),
            number = n())

sero_byage <- ggplot() + 
  geom_col(data = hsv_serobyage %>% 
               filter(taxon_species %in% c("Human herpesvirus 1", "Human herpesvirus 2",
                                           #"Varicella-zoster virus", 
                                           "Human herpesvirus 3",
                                           "Human cytomegalovirus", "Epstein-Barr virus",
                                           #"Human herpesvirus 6A", 
                                           "Human herpesvirus 6B")),
             aes(x = taxon_species, y = proportion,
                 fill = age_group),
             position = position_dodge(),
             color = "black"
             ) +
  ylab("Proportion positive") +
  scale_x_discrete(labels=c("Human herpesvirus 1" = "HSV1",
                            "Human herpesvirus 2" = "HSV2",
                            "Human herpesvirus 3" = "VZV",
                            "Epstein-Barr virus" = "EBV",
                            "Human cytomegalovirus" = "CMV",
                            "Human herpesvirus 6B" = "HHV6B"
  )) +
  theme_light() + 
  scale_fill_viridis_d() +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - child herpesvirus seropositivity by age.jpeg",
       sero_byage,
       scale = 1,
       height = 1.5,
       width = 3,
       units = "in",
       dpi = 300)

##############
## Cord vs child scores
##############
child_ages <- hsv_VAR %>% ungroup() %>% group_by(ID) %>%
  summarise(years2 = max(years))

VARscores_pairs <- hsv_VAR %>% ungroup() %>% left_join(child_ages, multiple = "all") %>%
  pivot_wider(names_from = Plasma_type, id_cols = c(taxon_species, taxon_genus, ID, years2),
              values_from = vir_score)

cord_v_child <- ggplot(data = VARscores_pairs %>% ungroup() %>% 
         filter(taxon_species %in% c("Human herpesvirus 1", "Human herpesvirus 2",
                                     #"Varicella-zoster virus", 
                                     "Human herpesvirus 3",
                                     "Human cytomegalovirus", "Epstein-Barr virus",
                                     #"Human herpesvirus 6A", 
                                     "Human herpesvirus 6B")) %>%
           mutate(`Age < 1 Year` = years2 < 1)) + 
  geom_point(aes(x = cord, y = child,
                 color = `Age < 1 Year`),
             size = .3,
             alpha = .5) +
  ylab("Child Reactivity") +
  xlab("Cord Reactivity") +
  theme_light() +
  scale_color_manual(values = c("black", "blue")) +
  #geom_smooth(aes(x = cord, y = child,
  #                color = years2 < 1),
  #            method = "lm",
  #            formula = y ~ x) +
  facet_wrap(vars(taxon_species)) + 
  theme(text = element_text(size = 8),
        legend.position = "none")

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - cord v child herpesvirus VARscore.jpeg",
       cord_v_child,
       scale = 1,
       height = 3,
       width = 4,
       units = "in",
       dpi = 300)

## remake for legend
cord_v_child <- ggplot(data = VARscores_pairs %>% ungroup() %>% 
                         filter(taxon_species %in% c("Human herpesvirus 1", "Human herpesvirus 2",
                                                     #"Varicella-zoster virus", 
                                                     "Human herpesvirus 3",
                                                     "Human cytomegalovirus", "Epstein-Barr virus",
                                                     #"Human herpesvirus 6A", 
                                                     "Human herpesvirus 6B")) %>%
                         mutate(`Age < 1 Year` = years2 < 1)) + 
  geom_point(aes(x = cord, y = child,
                 color = `Age < 1 Year`),
             size = .3,
             alpha = .5) +
  ylab("Child Reactivity") +
  xlab("Cord Reactivity") +
  theme_light() +
  scale_color_manual(values = c("black", "blue")) +
  #geom_smooth(aes(x = cord, y = child,
  #                color = years2 < 1),
  #            method = "lm",
  #            formula = y ~ x) +
  facet_wrap(vars(taxon_species)) + 
  theme(text = element_text(size = 8))

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - cord v child herpesvirus VARscore - legend.jpeg",
       cord_v_child,
       scale = 1,
       height = 3,
       width = 4,
       units = "in",
       dpi = 300)

##############
# cross reactivity 
##############

HSV1_HSV2 <- hsv_VAR %>% ungroup() %>%
  pivot_wider(id_cols = c(ID, Plasma_type), names_from = taxon_species,
              values_from = vir_score, values_fn = max)

line_1 <- data.frame(x = 2, y = 1, xend = 25, yend = 12.5)
line_2 <- data.frame(x = 1, y = 2, xend = 15, yend = 30)

HSV_xreact <- ggplot() + 
  geom_point(data = HSV1_HSV2,
             aes(x = `Human herpesvirus 1`, 
                 y = `Human herpesvirus 2`,
                 color = Plasma_type),
             size = .3,
             alpha = .3) +
  xlab("HSV1") +
  ylab("HSV2") +
  theme_light() +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(vars(Plasma_type)) +
    theme(text = element_text(size = 8),
          legend.position = "none") + 
  #geom_hline(yintercept = 1.936591)+ 
  #geom_vline(xintercept = 1.936591)+
  geom_segment(data = line_1,
               aes(x = x, y = y, xend = xend, yend = yend), color = "grey20")+
  geom_segment(data = line_2,
               aes(x = x, y = y, xend = xend, yend = yend), color = "grey20")

HSV1_HSV2 <- hsv_VAR %>% ungroup() %>%
  pivot_wider(id_cols = c(ID, Plasma_type), names_from = taxon_species,
              values_from = vir_score, values_fn = max)

HSV1_HSV2 <- HSV1_HSV2 %>% mutate(simplex = case_when(
  (`Human herpesvirus 1` < 1.936591) & (`Human herpesvirus 2` < 1.936591) ~ "Double Negative",
  (`Human herpesvirus 1` >= 1.936591) & ((`Human herpesvirus 2` < 1.936591) | (`Human herpesvirus 2` < 0.5*`Human herpesvirus 1`)) ~ "HSV1 Single Positive",
  (`Human herpesvirus 2` >= 1.936591) & ((`Human herpesvirus 1` < 1.936591) | (`Human herpesvirus 1` < 0.5*`Human herpesvirus 2`)) ~ "HSV2 Single Positive",
  T ~ "Double Positive"
)) %>% ungroup() %>%
  group_by(Plasma_type, simplex) %>% summarise(n = n())
  
HSV1_HSV2_sero <- HSV1_HSV2 %>% mutate(sero = n) %>% mutate(n = 965) %>%
  mutate(percent = sero / n)

HSV1_HSV2_sero_byage <- hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
  pivot_wider(id_cols = c(ID, years), names_from = taxon_species,
              values_from = vir_score, values_fn = max) %>% mutate(simplex = case_when(
  (`Human herpesvirus 1` < 1.936591) & (`Human herpesvirus 2` < 1.936591) ~ "Double Negative",
  (`Human herpesvirus 1` >= 1.936591) & ((`Human herpesvirus 2` < 1.936591) | (`Human herpesvirus 2` < 0.5*`Human herpesvirus 1`)) ~ "HSV1 Single Positive",
  (`Human herpesvirus 2` >= 1.936591) & ((`Human herpesvirus 1` < 1.936591) | (`Human herpesvirus 1` < 0.5*`Human herpesvirus 2`)) ~ "HSV2 Single Positive",
  T ~ "Double Positive"
)) %>% ungroup() %>% mutate(HSV1 = (simplex %in% c("HSV1 Single Positive", "Double Positive"))) %>%
  mutate(HSV2 = (simplex %in% c("HSV2 Single Positive", "Double Positive"))) %>%
  group_by(ID, years, simplex) %>%
mutate(age_group = case_when(
  years < 1 ~ "0 - 1",
  years < 2 ~ "1 - 2",
  years < 3 ~ "2 - 3",
  years < 4 ~ "3 - 4",
  years >= 4 ~ "4+"
)) %>% group_by(age_group) %>%
  summarise(`HSV1 proportion` = sum(HSV1)/n(), 
            `HSV1 positive` = sum(HSV1),
            number = n(),
            `HSV2 proportion` = sum(HSV2)/n(), 
            `HSV2 positive` = sum(HSV2))


ggsave("Figures/20240424 - HSV cross reactivity.jpeg",
       HSV_xreact,
       scale = 1,
       height = 1.5,
       width = 3,
       units = "in",
       dpi = 300)

HSV1_HSV2 <- hsv_VAR %>% ungroup() %>%
  pivot_wider(id_cols = c(ID, Plasma_type), names_from = taxon_species,
              values_from = vir_score, values_fn = max)

HSV6_xreact <- ggplot() + 
  geom_point(data = HSV1_HSV2,
             aes(x = `Human herpesvirus 6A`, 
                 y = `Human herpesvirus 6B`,
                 color = Plasma_type),
             size = .3,
             alpha = .3) +
  xlab("HHV 6a") +
  ylab("HHV 6b") +
  theme_light() +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(vars(Plasma_type)) +
  theme(text = element_text(size = 8),
        legend.position = "none") + 
  geom_hline(yintercept = 1.936591) + 
  geom_vline(xintercept = 1.936591) +
  geom_abline(intercept = 0, slope = 2)+
  geom_abline(intercept = 0, slope = 0.5)

HSV1_HSV2 <- HSV1_HSV2 %>% mutate(simplex = case_when(
  (`Human herpesvirus 6A` < 1.936591) & (`Human herpesvirus 6B` < 1.936591) ~ "Double Negative",
  (`Human herpesvirus 6A` >= 1.936591) & ((`Human herpesvirus 6B` < 1.936591) | (`Human herpesvirus 6B` < 0.5*`Human herpesvirus 6A`)) ~ "HSV1 Single Positive",
  (`Human herpesvirus 6B` >= 1.936591) & ((`Human herpesvirus 6A` < 1.936591) | (`Human herpesvirus 6A` < 0.5*`Human herpesvirus 6B`)) ~ "HSV2 Single Positive",
  T ~ "Double Positive"
)) %>% ungroup() %>%
  group_by(Plasma_type, simplex) %>% summarise(n = n())

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - HHV 6a v 6b cross reactivity.jpeg",
       HSV6_xreact,
       scale = 1,
       height = 1.5,
       width = 3,
       units = "in",
       dpi = 300)

HSV1_HSV2 <- HSV1_HSV2 %>% ungroup() %>%
  filter(Plasma_type == "child")

library(GGally)

HSV1_HSV2 <- hsv_VAR %>% ungroup() %>%
  pivot_wider(id_cols = c(ID, Plasma_type), names_from = taxon_species,
              values_from = vir_score, values_fn = max)

names(HSV1_HSV2) <- c("ID", "Plasma_type", 
                      "HSV1", "HSV2",
                      "VZV", "EBV",
                      "CMV", "HSV6A", "HSV6B",
                      "HSV7", "HSV8")

pairs_dots <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
    geom_point(..., alpha = 0.05) 
}

ggpairs(data = HSV1_HSV2,
        columns = c(3:11),
        #mapping = aes(alpha = .1)
        lower = list(continuous = pairs_dots)) +
  theme_light()
