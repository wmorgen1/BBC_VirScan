setwd("~/Larman_Lab")
library(tidyverse)
library(RColorBrewer)
library(gtable)
library(grid)
library(ggplot2)
library(ggbeeswarm)

library(data.table)
'%nin%' <- function(x,y)!('%in%'(x,y))

#library(ARscore)

########
# Data Prep
#######
library(readxl)
sample_meta <- read_excel("BBC/sample_meta.xls")

sample_meta <- sample_meta %>% mutate(age = difftime(time1 = Plasma_Date, time2 = DOD, units = "days")) %>%
  mutate(years = (age %>% as.numeric())/365.25) %>% mutate(Plasma_type = ifelse(Plasma_type=="Postnatal","child",Plasma_type))

breast_fed_meta <- read_excel("BBC/breastfed_update ID with codebook.xlsx")

sample_meta <- sample_meta %>% left_join(breast_fed_meta) %>% mutate(breast_fed1 = ifelse(breast_fed1 == "2", "1", breast_fed1))

BBC_VARscores <- read_rds("BBC/Data/20240424 - BBC VARscores - maternal abs removed.rds")

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
  pivot_wider(id_cols = c(ID, Plasma_type, breast_fed1), names_from = taxon_species,
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
  dplyr::select(`ID`, `Mom CMV positive`, `Mom HSV1 positive`, `Mom HSV2 positive`, `breast_fed1`)

dotplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
               left_join(Maternal_pos) %>% filter(taxon_species == "Human cytomegalovirus") %>%
               filter(`Mom CMV positive` == TRUE) %>% filter(breast_fed1 != "NA"),
             aes(x = `breast_fed1`, y = vir_score,
                 color = `breast_fed1`),
             position = position_quasirandom(),
             size = .1,
             alpha = .3) +
  ylab("Child CMV VARscore") +
  geom_hline(yintercept = 1.936591) +
  scale_x_discrete(labels=c(
    "1" = "Breastfed",
    "0" = "Formula"
    )) +
  theme_light() +
  scale_color_manual(values = c("black", "blue")) +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.position = "none")

ggsave("BBC/Figures/20250329 - CMV VARscore dotplot by breastfeeding.jpeg",
       dotplot,
       scale = 1,
       height = 1.5,
       width = 1.5,
       units = "in",
       dpi = 300)

ftest_cmv <- hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
  left_join(Maternal_pos) %>% filter(taxon_species == "Human cytomegalovirus") %>%
  filter(`Mom CMV positive` == TRUE) %>% filter(breast_fed1 != "NA") %>%
  mutate(`kid pos` = vir_score > 1.936591) %>%
  group_by(`breast_fed1`, `kid pos`) %>% summarise(n = n())

#### old groups
# 163- , 21+  formula
# 361- , 229+ mixed
# 30-  , 31+  breast
#### New groups
# 163- , 21+  formula
# 391- , 260 breast


fisher.test(matrix(data = c(260, 391,
                            21,  163), nrow = 2)) # any breastfeeding vs formula only

timeplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
               left_join(Maternal_pos) %>% filter(taxon_species == "Human cytomegalovirus") %>%
               filter(`Mom CMV positive` == TRUE) %>% filter(breast_fed1 != "NA") %>%
               mutate(breast_fed1 = case_when(
                 breast_fed1 == "0" ~ "Formula",
                 breast_fed1 == "1" ~ "Breastfed"
                 )) %>% mutate(breast_fed1 = as.factor(breast_fed1)) %>%
               mutate(breast_fed1 = fct_relevel(breast_fed1, c("Formula", "Breastfed"))),
             aes(x = years, y = vir_score,
                 color = `breast_fed1`),
             size = .3,
             alpha = .3) +
  geom_hline(yintercept = 1.936591) +
  ylab("VARscore") +
  xlab("Age (years)") +
  theme_light() + 
  scale_color_manual(values = c("black", "blue")) +
  facet_wrap(vars(`breast_fed1`)) + 
  theme(text = element_text(size = 8),
        legend.position = "none")

ggsave("BBC/Figures/20250329 - child CMV VARscore vs time by feeding dotplot.jpeg",
       timeplot,
       scale = 1,
       height = 1.5,
       width = 3,
       units = "in",
       dpi = 300)

#### HSV 1

dotplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
               left_join(Maternal_pos) %>% filter(taxon_species == "Human herpesvirus 1") %>%
               filter(`Mom HSV1 positive` == TRUE) %>% filter(breast_fed1 != "NA"),
             aes(x = `breast_fed1`, y = vir_score,
                 color = `breast_fed1`),
             position = position_quasirandom(),
             size = .1,
             alpha = .3) +
  ylab("Child HSV1 VARscore") +
  geom_hline(yintercept = 1.936591) +
  scale_x_discrete(labels=c(
    "1" = "Breastfed",
    "0" = "Formula"
  )) +
  theme_light() +
  scale_color_manual(values = c("black", "blue")) +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank(),
        legend.position = "none")

ggsave("BBC/Figures/20250329 - HSV1 VARscore dotplot by breastfeeding.jpeg",
       dotplot,
       scale = 1,
       height = 1.5,
       width = 1.5,
       units = "in",
       dpi = 300)

ftest_hsv <- hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
  left_join(Maternal_pos) %>% filter(taxon_species == "Human herpesvirus 1") %>%
  filter(`Mom HSV1 positive` == TRUE) %>% filter(breast_fed1 != "NA") %>%
  mutate(`kid pos` = vir_score > 1.936591) %>%
  group_by(`breast_fed1`, `kid pos`) %>% summarise(n = n())

#### old groups
# 167- , 23+  formula
# 493- , 74+ mixed
# 47-  , 6+  breast
#### new groups
# 167- , 23+  formula
# 540- , 80+  breast

fisher.test(matrix(data = c(167, 23,
                            540,  80), nrow = 2)) # any breastfeeding vs formula only

timeplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
               left_join(Maternal_pos) %>% filter(taxon_species == "Human herpesvirus 1") %>%
               filter(`Mom HSV1 positive` == TRUE) %>% filter(breast_fed1 != "NA") %>%
               mutate(breast_fed1 = case_when(
                 breast_fed1 == "0" ~ "Formula",
                 breast_fed1 == "1" ~ "Breastfed"
               )) %>% mutate(breast_fed1 = as.factor(breast_fed1)) %>%
               mutate(breast_fed1 = fct_relevel(breast_fed1, c("Formula", "Breastfed"))),
             aes(x = years, y = vir_score,
                 color = `breast_fed1`),
             size = .3,
             alpha = .3) +
  geom_hline(yintercept = 1.936591) +
  ylab("VARscore") +
  xlab("Age (years)") +
  theme_light() + 
  scale_color_manual(values = c("black", "blue")) +
  facet_wrap(vars(`breast_fed1`)) + 
  theme(text = element_text(size = 8),
        legend.position = "none")

ggsave("BBC/Figures/20250329 - child HSV1 VARscore vs time by feeding dotplot.jpeg",
       timeplot,
       scale = 1,
       height = 1.5,
       width = 3,
       units = "in",
       dpi = 300)

