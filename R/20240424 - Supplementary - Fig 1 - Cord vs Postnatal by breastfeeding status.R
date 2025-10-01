setwd("C:/Users/wmorgen1/Documents/Larman_Lab/")
library(tidyverse)
library(lubridate)
library(RColorBrewer)
library(gtable)
library(grid)
library(ggplot2)
library(ggbeeswarm)
library(data.table)
#library(ARscore)

#library(GGally)
library(gtable)
library(grid)

########
# Data Prep
#######
library(readxl)
sample_meta <- read_excel("BBC/sample_meta.xls")
#View(sample_meta)
sample_meta <- sample_meta %>% mutate(age = difftime(time1 = Plasma_Date, time2 = DOD, units = "days")) %>%
  mutate(years = (age %>% as.numeric())/365.25)

sample_meta <- sample_meta %>% filter(lab_id %in% names(BBC_hfc))

IDs_torun <- sample_meta %>% group_by(ID) %>% summarise(n = n()) %>%
  filter(n == 2)

phiphfc_1 <- fread("C:/Users/wmorgen1/Documents/Larman_Lab/BBC/Data/phipseq_0232_0241_VirscanLar_000_Hits_foldchange_annotated.tsv") %>%
  filter(taxon_species %in% BBC_VARscores_pairs_mat_removed_totest$taxon_species)
phiphfc_2 <- fread("C:/Users/wmorgen1/Documents/Larman_Lab/BBC/Data/phipseq_0242_0254_VirscanLar_000_Hits_foldchange_annotated.tsv") %>%
  filter(taxon_species %in% BBC_VARscores_pairs_mat_removed_totest$taxon_species)

phiphfc <- phiphfc_1 %>% full_join(phiphfc_2)

BBC_hfc <- phiphfc %>% dplyr::select(-contains("Frank")) %>% 
  dplyr::select(-contains("control")) %>% 
  dplyr::select(-contains("Control")) %>% 
  dplyr::select(-contains("Beads")) %>% 
  dplyr::select(-contains("BEADS")) %>% 
  dplyr::select(-contains("beads"))

names(BBC_hfc) <- names(BBC_hfc) %>% str_replace_all("-", ".")
names(BBC_hfc) <- paste0("X", names(BBC_hfc))

library(readxl)
sample_meta <- read_excel("BBC/sample_meta.xls")
#View(sample_meta)
sample_meta <- sample_meta %>% mutate(age = difftime(time1 = Plasma_Date, time2 = DOD, units = "days")) %>%
  mutate(years = (age %>% as.numeric())/365.25)

correlation_data <- data.frame(matrix(data = 0, nrow = length(sample_meta$ID %>% unique()), ncol = 5))
names(correlation_data) <- c("ID",
                             "Child Age",
                             "Pearson R",
                             "Spearman rho",
                             "peptide lm coefficient")
IDs <- (sample_meta$ID %>% unique())

# not 30, 35, 49, 53
for(i in c(1:965)){
  correlation_data[i, 1] <- IDs[i]
  correlation_data[i, 2] <- (sample_meta %>% filter(Plasma_type == "Postnatal") %>%
                               filter(ID == IDs[i]))$years
  
  curr_correlation <- cor(x = (BBC_hfc %>% select(sample_meta %>% filter(ID == IDs[i]) %>%
                                                    filter(Plasma_type == "Postnatal") %>%
                                                    dplyr::select(lab_id) %>%
                                                    pull)),
                          y = (BBC_hfc %>% select(sample_meta %>% filter(ID == IDs[i]) %>%
                                                    filter(Plasma_type == "cord") %>%
                                                    dplyr::select(lab_id) %>%
                                                    pull)),
                          method = "pearson")
  
  correlation_data[i, 3] <- curr_correlation
  
  curr_correlation <- cor(x = (BBC_hfc %>% select(sample_meta %>% filter(ID == IDs[i]) %>%
                                                    filter(Plasma_type == "Postnatal") %>%
                                                    dplyr::select(lab_id) %>%
                                                    pull)),
                          y = (BBC_hfc %>% select(sample_meta %>% filter(ID == IDs[i]) %>%
                                                    filter(Plasma_type == "cord") %>%
                                                    dplyr::select(lab_id) %>%
                                                    pull)),
                          method = "spearman")
  
  correlation_data[i, 4] <- curr_correlation
  
  curr_correlation <- lm(log2((BBC_hfc %>% select(sample_meta %>% filter(ID == IDs[i]) %>%
                                                    filter(Plasma_type == "Postnatal") %>%
                                                    dplyr::select(lab_id) %>%
                                                    pull) %>% pull)) ~
                           log2((BBC_hfc %>% select(sample_meta %>% filter(ID == IDs[i]) %>%
                                                      filter(Plasma_type == "cord") %>%
                                                      dplyr::select(lab_id) %>%
                                                      pull) %>% pull)))
  
  correlation_data[i, 5] <- curr_correlation$coefficients[2]
  print(i)
}

#write.csv(correlation_data, "/BBC/20240315 - BBC cord child correlation.csv")
correlation_data <- read.csv("C:/Users/wmorgen1/Documents/Larman_Lab/BBC/Data/20240315 - BBC cord child correlation.csv")

breast_fed_meta <- read_excel("BBC/breastfed_update ID with codebook.xlsx")

names(correlation_data) <- names(correlation_data) %>% str_replace_all("\\.", " ")
breastfed <- 
#pearson by age
corr_plot = ggplot(data = correlation_data %>% left_join(breast_fed_meta),
                   aes(y = `Pearson R`,
                       x = `Child Age`,
                       color = `breast_fed1`)) +
  geom_point(alpha = .5, size = .5) +
  theme_light() +
  #geom_smooth() +
  theme(text = element_text(size = 8))

a <- correlation_data %>% left_join(breast_fed_meta) %>% filter(`Child Age` < 2)

corr_plot = ggplot(data = a,
                   aes(y = `Pearson R`,
                       x = `Child Age`,
                       color = `breast_fed1`)) +
  geom_point(alpha = .5, size = 2.5) +
  theme_light() +
  geom_smooth() +
  theme(text = element_text(size = 8))


#ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240315 - cord child pearson vs age.jpeg",
#       corr_plot,
#       scale = 1,
#       height = 2,
#       width = 2,
#       units = "in",
#       dpi = 300) 