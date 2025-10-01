setwd("~/R/Larman_Lab")
library(tidyverse)
library(lubridate)
library(RColorBrewer)
library(gtable)
library(grid)
library(ggplot2)
library(ggbeeswarm)
source("ArboScan/common_functions.R")
library(data.table)
library(viridis)
#library(ARscore)

########
# Data Prep
#######
# row 1153 id 5272 has age wrong. I edited to <1 yo instead of -6 could exclude
phip_VARscores_1 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 1.rds")
phip_VARscores_2 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 2.rds")
phip_VARscores_3 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 3.rds")
phip_VARscores_4 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 4.rds")
phip_VARscores_5 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 5.rds")
phip_VARscores_6 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 6.rds")
phip_VARscores_7 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 7.rds")
phip_VARscores_8 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 8.rds")
phip_VARscores_9 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 9.rds")
phip_VARscores_10 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 10.rds")
phip_VARscores_11 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 11.rds")

BBC_VARscores <- phip_VARscores_1 %>% full_join(phip_VARscores_2) %>% full_join(phip_VARscores_3) %>% 
  full_join(phip_VARscores_4) %>% full_join(phip_VARscores_5) %>% full_join(phip_VARscores_6) %>% 
  full_join(phip_VARscores_7) %>% full_join(phip_VARscores_8) %>% full_join(phip_VARscores_9) %>% 
  full_join(phip_VARscores_10) %>% full_join(phip_VARscores_11) 

BBC_VARscores <- BBC_VARscores %>% filter(!grepl("Frank", sample_id)) %>% filter(!grepl("control", sample_id)) %>%
  filter(!grepl("Control", sample_id))

BBC_VARscores <- BBC_VARscores %>% mutate(sample_id = str_replace_all(sample_id, "-", ".")) %>%
  mutate(sample_id = paste0("X", sample_id))

library(readxl)
sample_meta <- read_excel("BBC/sample_meta.xls")
View(sample_meta)
sample_meta <- sample_meta %>% mutate(age = difftime(time1 = Plasma_Date, time2 = DOD, units = "days")) %>%
  mutate(years = (age %>% as.numeric())/365.25)

BBC_VARscores <- BBC_VARscores %>% filter(sample_id %in% sample_meta$lab_id) %>%
  rename(lab_id = sample_id) %>%
  left_join(sample_meta, multiple = "all")

# Set positives and negatives

pos_control_viruses = c("Enterovirus B", "Epstein-Barr virus", 
                        "Human respiratory syncytial virus A")
neg_control_viruses = c("Ebolavirus", "Coltivirus", "Marburgvirus",
                        "Thogotovirus", "Hepacivirus", "Yatapoxvirus",
                        "Lentivirus")

library(pROC)
posneg <- BBC_VARscores %>% ungroup() %>% filter(Plasma_type == "cord") %>%
  mutate(positive = ifelse(taxon_species %in% pos_control_viruses, "Pos", "Other")) %>%
  mutate(positive = ifelse(taxon_genus %in% neg_control_viruses, "Neg", positive)) %>%
  #left_join(number_hits) %>% mutate(positive = ifelse(hits < 3, "Neg", positive)) %>%
  filter(positive %in% c("Neg", "Pos")) 

z <- roc(response = posneg$positive, predictor = posneg$vir_score, levels = c("Pos", "Neg"))
z
plot(z)
coords(z, x = "best", best.method = 'y')

ggplot() + 
  geom_point(data = posneg %>% ungroup() %>% group_by(Plasma_type),
             aes(x = taxon_species, y = vir_score,
                 color = Plasma_type),
             position = position_quasirandom(),
             size = 1,
             alpha = .3) +
  ylab("VARscore") +
  xlab("Virus") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  geom_hline(yintercept = 1.936591) 

# check means and standard deviations for each virus
BBC_VARscores_meansd <- BBC_VARscores %>% filter(Plasma_type == "cord") %>%
  filter(vir_score < (1.936591)) %>%
  group_by(taxon_genus, taxon_species) %>% summarise(mean_score = mean(vir_score),
                                                     score_sd = sd(vir_score)) %>%
  mutate(threshold = mean_score + 3 * score_sd)

#1.936591