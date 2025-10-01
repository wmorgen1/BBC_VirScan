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
  filter(taxon_species %in% c("Human herpesvirus 1", "Human herpesvirus 2",
                              #"Varicella-zoster virus", 
                              "Human herpesvirus 3",
                              "Human cytomegalovirus", "Epstein-Barr virus",
                              "Human herpesvirus 6A", "Human herpesvirus 6B",
                              "Human herpesvirus 7", "Human herpesvirus 8")) %>%
  mutate(taxon_species = as.factor(taxon_species)) %>%
  mutate(taxon_species = fct_relevel(taxon_species, "Human herpesvirus 1", "Human herpesvirus 2",
                                     #"Varicella-zoster virus", 
                                     "Human herpesvirus 3",
                                     "Epstein-Barr virus", "Human cytomegalovirus",
                                     "Human herpesvirus 6A", "Human herpesvirus 6B",
                                     "Human herpesvirus 7", "Human herpesvirus 8")) %>%
  arrange(taxon_species)

Maternal_pos <- hsv_VAR %>% ungroup() %>%
  pivot_wider(id_cols = c(ID, Plasma_type), names_from = taxon_species,
              values_from = vir_score, values_fn = max) %>%
  filter(Plasma_type == "cord") %>%
  mutate(`Mom CMV positive` = `Human cytomegalovirus` > 1.936591) %>%
  mutate(simplex = case_when(
    (`Human herpesvirus 1` < 1.936591) & (`Human herpesvirus 2` < 1.936591) ~ "Double Negative",
    (`Human herpesvirus 1` >= 1.936591) & ((`Human herpesvirus 2` < 1.936591) | (`Human herpesvirus 2` < 0.5*`Human herpesvirus 1`)) ~ "HSV1 Single Positive",
    (`Human herpesvirus 2` >= 1.936591) & ((`Human herpesvirus 1` < 1.936591) | (`Human herpesvirus 1` < 0.5*`Human herpesvirus 2`)) ~ "HSV2 Single Positive",
    T ~ "Double Positive"
  )) %>% mutate(`Mom HSV1 positive` = simplex %in% c("Double Positive", "HSV1 Single Positive")) %>%
  mutate(`Mom HSV2 positive` = simplex %in% c("Double Positive", "HSV2 Single Positive")) %>%
  dplyr::select(`ID`, `Mom CMV positive`, `Mom HSV1 positive`, `Mom HSV2 positive`)

dotplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
               left_join(Maternal_pos) %>% filter(taxon_species == "Human cytomegalovirus") %>%
               filter(!is.na(`Mom CMV positive`)),
             aes(x = `Mom CMV positive`, y = vir_score,
                 color = `Mom CMV positive`),
             position = position_quasirandom(),
             size = .1,
             alpha = .3) +
  ylab("Child CMV VARscore") +
  geom_hline(yintercept = 1.936591) +
  scale_x_discrete(labels=c("TRUE" = "Mom\nCMV+",
                            "FALSE" = "Mom\nCMV-"
                            )) +
  theme_light() +
  scale_color_manual(values = c("black", "darkgreen")) +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.position = "none")

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240515 - CMV VARscore dotplot by maternal status.jpeg",
       dotplot,
       scale = 1,
       height = 1.5,
       width = 1.5,
       units = "in",
       dpi = 300)

dotplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
               left_join(Maternal_pos) %>% filter(taxon_species == "Human herpesvirus 1") %>%
               filter(!is.na(`Mom HSV1 positive`)),
             aes(x = `Mom HSV1 positive`, y = vir_score,
                 color = `Mom HSV1 positive`),
             position = position_quasirandom(),
             size = .1,
             alpha = .3) +
  ylab("Child HSV1 VARscore") +
  geom_hline(yintercept = 1.936591) +
  scale_x_discrete(labels=c("TRUE" = "Mom\nHSV1+",
                            "FALSE" = "Mom\nHSV1-"
  )) +
  theme_light() +
  scale_color_manual(values = c("black", "darkgreen")) +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.position = "none")

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240515 - HSV1 VARscore dotplot by maternal status.jpeg",
       dotplot,
       scale = 1,
       height = 1.5,
       width = 1.5,
       units = "in",
       dpi = 300)

ftest_HSV1 <- hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
  left_join(Maternal_pos) %>% filter(taxon_species == "Human herpesvirus 1") %>%
  filter(!is.na(`Mom CMV positive`)) %>% mutate(`kid pos` = vir_score > 1.936591) %>%
  group_by(`Mom HSV1 positive`, `kid pos`) %>% summarise(n = n())

fisher.test(matrix(data = c(104, 720,
                            6, 135), nrow = 2))

ftest_cmv <- hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
  left_join(Maternal_pos) %>% filter(taxon_species == "Human cytomegalovirus") %>%
  filter(!is.na(`Mom CMV positive`)) %>% mutate(`kid pos` = vir_score > 1.936591) %>%
  group_by(`Mom CMV positive`, `kid pos`) %>% summarise(n = n())

fisher.test(matrix(data = c(285, 563,
                            13, 104), nrow = 2))

timeplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
               left_join(Maternal_pos) %>% filter(taxon_species == "Human cytomegalovirus") %>%
               filter(!is.na(`Mom CMV positive`)) %>%
               mutate(`Mom CMV positive` = ifelse(`Mom CMV positive`, "Mom CMV+", "Mom CMV-")),
             aes(x = years, y = vir_score,
                 color = `Mom CMV positive`),
             size = .3,
             alpha = .3) +
  geom_hline(yintercept = 1.936591) +
  ylab("VARscore") +
  xlab("Age (years)") +
  theme_light() + 
  scale_color_manual(values = c("black", "darkgreen")) +
  facet_wrap(vars(`Mom CMV positive`)) + 
  theme(text = element_text(size = 8),
        legend.position = "none")

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240515 - child CMV VARscore vs time by mom serostatus dotplot.jpeg",
       timeplot,
       scale = 1,
       height = 1.5,
       width = 3,
       units = "in",
       dpi = 300)

timeplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
               left_join(Maternal_pos) %>% filter(taxon_species == "Human herpesvirus 1") %>%
               filter(!is.na(`Mom HSV1 positive`)) %>%
               mutate(`Mom HSV1 positive` = ifelse(`Mom HSV1 positive`, "Mom HSV1+", "Mom HSV1-")),
             aes(x = years, y = vir_score,
                 color = `Mom HSV1 positive`),
             size = .3,
             alpha = .3) +
  geom_hline(yintercept = 1.936591) +
  ylab("VARscore") +
  xlab("Age (years)") +
  theme_light() + 
  scale_color_manual(values = c("black", "darkgreen")) +
  facet_wrap(vars(`Mom HSV1 positive`)) + 
  theme(text = element_text(size = 8),
        legend.position = "none")

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240515 - child HSV1 VARscore vs time by mom serostatus dotplot.jpeg",
       timeplot,
       scale = 1,
       height = 1.5,
       width = 3,
       units = "in",
       dpi = 300)
