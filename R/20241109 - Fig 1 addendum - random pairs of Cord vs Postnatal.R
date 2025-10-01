setwd("C:/Users/wmorgen1/Documents/Larman_Lab/BBC/")
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

viruses_totest <- read.csv("Data/20240424 - BBC cord child volcano.csv")

phiphfc_1 <- fread("Data/phipseq_0232_0241_VirscanLar_000_Hits_foldchange_annotated.tsv") %>%
  filter(taxon_species %in% viruses_totest$virus)
phiphfc_2 <- fread("Data/phipseq_0242_0254_VirscanLar_000_Hits_foldchange_annotated.tsv") %>%
  filter(taxon_species %in% viruses_totest$virus)

phiphfc <- phiphfc_1 %>% full_join(phiphfc_2)

BBC_hfc <- phiphfc %>% dplyr::select(-contains("Frank")) %>% 
  dplyr::select(-contains("control")) %>% 
  dplyr::select(-contains("Control")) %>% 
  dplyr::select(-contains("Beads")) %>% 
  dplyr::select(-contains("BEADS")) %>% 
  dplyr::select(-contains("beads")) %>% 
  dplyr::select(-contains("Empty"))

names(BBC_hfc) <- names(BBC_hfc) %>% str_replace_all("-", ".")
names(BBC_hfc) <- paste0("X", names(BBC_hfc))

library(readxl)
sample_meta <- read_excel("sample_meta.xls")
#View(sample_meta)
sample_meta <- sample_meta %>% mutate(age = difftime(time1 = Plasma_Date, time2 = DOD, units = "days")) %>%
  mutate(years = (age %>% as.numeric())/365.25)

sample_meta <- sample_meta %>% filter(lab_id %in% names(BBC_hfc))

IDs_torun <- sample_meta %>% group_by(ID) %>% summarise(n = n()) %>%
  filter(n == 2)

sample_meta <- sample_meta %>% filter(ID %in% IDs_torun$ID)

BBC_hfc <- BBC_hfc %>% dplyr::select(sample_meta$lab_id)

correlation_data <- cor(BBC_hfc, method = "spearman")

write.csv(correlation_data, "Data/20250330 - peptide peptide correlations all samples.csv")

names(correlation_data) <- row.names(correlation_data)
cord <- sample_meta %>% filter(Plasma_type == "cord") %>%
  dplyr::select(lab_id, ID) %>% dplyr::rename(cord = lab_id, cord_ID = ID)
child <- sample_meta %>% filter(Plasma_type == "Postnatal") %>%
  dplyr::select(lab_id, ID, years)

cor_data_to_plot <- correlation_data %>% as.data.frame() %>% 
  dplyr::select(cord$cord) %>%
  rownames_to_column() %>%
  dplyr::rename(lab_id = rowname) %>%
  filter(lab_id %in% child$lab_id) %>%
  pivot_longer(cols = -lab_id, names_to = "cord") %>%
  left_join(cord) %>% left_join(child) %>% filter(ID != cord_ID) %>%
  dplyr::rename(`Spearman rho` = value, `Child Age` = years) %>%
  group_by(`Child Age`, `lab_id`) %>% 
  summarise(`Median Spearman rho` = median(`Spearman rho`),
            `Mean Spearman rho` = mean(`Spearman rho`))

#spearman by age
corr_plot = ggplot(data = cor_data_to_plot,
                   aes(y = `Median Spearman rho`,
                       x = `Child Age`)) +
  geom_point(alpha = .5, size = .5) +
  theme_light() +
  geom_smooth() +
  theme(text = element_text(size = 8))

ggsave("Figures/20250330 - median nonpaired cord child spearman vs age.jpeg",
       corr_plot,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300) 

correlation_data_pair <- read.csv("Data/20240315 - BBC cord child correlation.csv")

correlation_data_pair <- correlation_data_pair %>% dplyr::select(-X, -Pearson.R, -V5) %>%
  left_join(sample_meta) %>% filter(Plasma_type == "Postnatal") %>%
  dplyr::select(ID, Spearman.rho, lab_id) %>% left_join(cor_data_to_plot) %>%
  pivot_longer(cols = c("Spearman.rho", "Median Spearman rho")) %>%
  mutate(name = case_when(
    name == "Spearman.rho" ~ "Paired",
    name == "Median Spearman rho" ~ "Median of\nNonpaired"
  ))

#pearson by age
dot_plot_comp_corr = ggplot(data = correlation_data_pair %>%
                              filter(`Child Age` > 2),
                   aes(y = value,
                       x = name)) +
  geom_boxplot(alpha = .5, size = .5,
               linewidth = 1) +
  geom_path(aes(group = ID),
            linewidth = 0.3, alpha = 0.2)+
  theme_light() +
  theme(text = element_text(size = 8))

ggsave("Figures/20250330 - cor with nonpaired vs paired cord blood.jpeg",
       dot_plot_comp_corr,
       scale = 1,
       height = 2,
       width = 1.5,
       units = "in",
       dpi = 300) 

correlation_data_pair <- read.csv("Data/20240315 - BBC cord child correlation.csv")

correlation_data_pair <- correlation_data_pair %>% dplyr::select(-X, -Pearson.R, -V5) %>%
  left_join(sample_meta) %>% filter(Plasma_type == "Postnatal") %>%
  dplyr::select(ID, Spearman.rho, lab_id) %>% left_join(cor_data_to_plot)

t.test(correlation_data_pair$Spearman.rho, correlation_data_pair$`Median Spearman rho`,
       paired = T)