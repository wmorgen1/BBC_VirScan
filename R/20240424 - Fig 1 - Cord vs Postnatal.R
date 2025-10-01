setwd("~/Larman_Lab")
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
'%nin%' <- function(x,y)!('%in%'(x,y))

########
# Data Prep
#######
phip_VARscores_1 <- read_rds("BBC/Data/20240206 - VARscores - 1.rds")
phip_VARscores_2 <- read_rds("BBC/Data/20240206 - VARscores - 2.rds")
phip_VARscores_3 <- read_rds("BBC/Data/20240206 - VARscores - 3.rds")
phip_VARscores_4 <- read_rds("BBC/Data/20240206 - VARscores - 4.rds")
phip_VARscores_5 <- read_rds("BBC/Data/20240206 - VARscores - 5.rds")
phip_VARscores_6 <- read_rds("BBC/Data/20240206 - VARscores - 6.rds")
phip_VARscores_7 <- read_rds("BBC/Data/20240206 - VARscores - 7.rds")
phip_VARscores_8 <- read_rds("BBC/Data/20240206 - VARscores - 8.rds")
phip_VARscores_9 <- read_rds("BBC/Data/20240206 - VARscores - 9.rds")
phip_VARscores_10 <- read_rds("BBC/Data/20240206 - VARscores - 10.rds")
phip_VARscores_11 <- read_rds("BBC/Data/20240206 - VARscores - 11.rds")

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
#View(sample_meta)
sample_meta <- sample_meta %>% mutate(age = difftime(time1 = Plasma_Date, time2 = DOD, units = "days")) %>%
  mutate(years = (age %>% as.numeric())/365.25)

BBC_VARscores <- BBC_VARscores %>% filter(sample_id %in% sample_meta$lab_id) %>%
  rename(lab_id = sample_id) %>%
  left_join(sample_meta, multiple = "all")

##############
# Determine how much maternal reactivity to subtract from child
##############

BBC_VARscores_pairs <- BBC_VARscores %>% ungroup() %>%
  pivot_wider(names_from = Plasma_type, id_cols = c(taxon_species, taxon_genus, ID),
              values_from = vir_score)

BBC_VARscores_pairs <- BBC_VARscores_pairs %>% filter(!is.na(Postnatal)) %>%
  filter(!is.na(cord))

correlation_data <- data.frame(matrix(data = 0, nrow = length(BBC_VARscores_pairs$ID %>% unique()), ncol = 5))
names(correlation_data) <- c("ID",
                             "Child Age",
                             "Pearson R",
                             "Spearman rho",
                             "VARscore lm coefficient")
IDs <- (BBC_VARscores_pairs$ID %>% unique())

# not 30, 35, 49, 53
for(i in c(1:965)){
  correlation_data[i, 1] <- IDs[i]
  correlation_data[i, 2] <- (sample_meta %>% filter(Plasma_type == "Postnatal") %>%
                               filter(ID == IDs[i]))$years
  
  curr_correlation <- cor(x = (BBC_VARscores_pairs %>% filter(ID == IDs[i]) %>%
                                                    dplyr::select(Postnatal) %>%
                                                    pull),
                          y = (BBC_VARscores_pairs %>% filter(ID == IDs[i]) %>%
                                                    dplyr::select(cord) %>%
                                                    pull),
                          method = "pearson")
  
  correlation_data[i, 3] <- curr_correlation
  
  curr_correlation <- cor(x = (BBC_VARscores_pairs %>% filter(ID == IDs[i]) %>%
                                                    dplyr::select(Postnatal) %>%
                                                    pull),
                          y = (BBC_VARscores_pairs %>% filter(ID == IDs[i]) %>%
                                                    dplyr::select(cord) %>%
                                                    pull),
                          method = "spearman")
  
  correlation_data[i, 4] <- curr_correlation
  
  curr_correlation <- lm((BBC_VARscores_pairs %>% filter(ID == IDs[i]) %>%
                                                    dplyr::select(Postnatal) %>%
                                                    pull) ~
                         0 + (BBC_VARscores_pairs %>% filter(ID == IDs[i]) %>%
                                                      dplyr::select(cord) %>%
                                                      pull))
  
  correlation_data[i, 5] <- curr_correlation$coefficients[1]
  print(i)
}


#write.csv(correlation_data, "BBC/Data/20240315 - BBC cord child correlation - VARscore.csv")
correlation_data <- read.csv("BBC/Data/20240315 - BBC cord child correlation - VARscore.csv")

reg_plot = ggplot(data = correlation_data,
                  aes(y = `Pearson R`,
                      x = `Child Age`)) + 
  geom_point(alpha = .5, size = .5) +
  theme_light() +
  geom_smooth()

reg_plot = ggplot(data = correlation_data,
                  aes(y = `VARscore lm coefficient`,
                      x = `Child Age`)) + 
  geom_point(alpha = .5, size = .5) +
  theme_light() +
  geom_smooth()

### modeling the regression coefficient
mat_ab_model <- nls(`VARscore lm coefficient` ~ asym * (exp(`Child Age` * rate2) - 1) / (exp(`Child Age` * rate2) + 1),
                    data = correlation_data,
                    start = list(asym = 1,
                                 rate2 = 1.442695))

mat_ab_model <- nls(`VARscore lm coefficient` ~ exp(-`Child Age` / rate1) + 
                      asym * (exp(`Child Age` / rate2) - 1) / (exp(`Child Age` / rate2) + 1),
                    data = correlation_data,
                    start = list(asym = 1,
                                 rate1 = .1442695 * 2, # 1 month
                                 rate2 = 1.442695), # 1 year
                    )

#plot(mat_ab_model)
mat_ab_model
x <- seq(from = 0, to = 5, by = .1) %>% as.data.frame()
names(x) = "Child Age"
y <- predict(mat_ab_model, x)
plot(x$`Child Age`, y)
model_data <- x %>% mutate(model = y) %>% mutate(decay = exp(-`Child Age`/0.2257)) %>%
  mutate(growth = 0.9609 * (exp(`Child Age`/0.7660) - 1)/(exp(`Child Age`/0.7660) + 1))

reg_plot = ggplot(data = correlation_data,
                  aes(y = `VARscore lm coefficient`,
                      x = `Child Age`)) + 
  geom_point(alpha = .5, size = .5) +
  theme_light() +
  geom_smooth() +
  geom_path(data = model_data,
            aes(y = `model`,
                x = `Child Age`),
            color = "red")+
  geom_path(data = model_data,
            aes(y = `decay`,
                x = `Child Age`),
            color = "red")+
  geom_path(data = model_data,
            aes(y = `growth`,
                x = `Child Age`),
            color = "red") + 
  ylab("Linear Model Coefficient")


ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - cord child lm coef vs age.jpeg",
       reg_plot,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300)

#decay half life of 1.87732 months
#growth half life of 6.37141 months


#######
# overal seroprevelance and largest differences between cord and child
#######
#965 pairs of cord and child

child_age <- sample_meta %>% filter(Plasma_type == "Postnatal") %>%
  dplyr::select(ID, years)

BBC_VARscores_pairs <- BBC_VARscores %>% ungroup() %>%
  pivot_wider(names_from = Plasma_type, id_cols = c(taxon_species, taxon_genus, ID),
              values_from = vir_score)

BBC_VARscores_pairs <- BBC_VARscores_pairs %>% filter(!is.na(Postnatal)) %>%
  filter(!is.na(cord))

BBC_VARscores_pairs <- BBC_VARscores_pairs %>% left_join(child_age)

BBC_VARscores_pairs_mat_removed <- BBC_VARscores_pairs %>% mutate(mat_coef = exp(-years/0.2257)) %>%
  mutate(child = Postnatal - cord * mat_coef)
  
BBC_VARscores_mat_removed <- BBC_VARscores_pairs_mat_removed %>% dplyr::select(taxon_species, taxon_genus, ID, cord, child) %>%
  pivot_longer(cols = c(cord, child), names_to = "Plasma_type")

write_rds(BBC_VARscores_mat_removed, "/BBC/20240424 - BBC VARscores - maternal abs removed.rds")

#####

# only looking at viruses with > 5% positivity in whole cohort
BBC_VARscores_pairs_mat_removed_totest <- BBC_VARscores %>% ungroup() %>% 
  group_by(taxon_species) %>% mutate(pos = vir_score > 1.936591) %>%
  summarise(n = n(), positive = sum(pos)) %>%
  filter(positive >= 97)

volcano_data <- data.frame(matrix(data = 0, nrow = 68, ncol = 14))
names(volcano_data) <- c("virus",
                         "Spearman rho",
                         "Spearman p",
                         "Pearson R",
                         "Pearson p",
                         "Wilcox p",
                         "ttest p",
                         "Difference in means",
                         "cord pos",
                         "cord neg",
                         "child pos",
                         "child neg",
                         "odds ratio",
                         "fisher p")
viruses <- (BBC_VARscores_pairs_mat_removed_totest$taxon_species %>% unique())

for(i in c(1:68)){
  volcano_data[i, 1] <- viruses[i]
  
  curr_spearman <- cor.test((BBC_VARscores_pairs_mat_removed_totest %>% ungroup() %>%
                               filter(taxon_species == viruses[i]))$cord,
                            (BBC_VARscores_pairs_mat_removed_totest %>% ungroup() %>%
                               filter(taxon_species == viruses[i]))$child,
                            method = "spearman")
  
  volcano_data[i, 2] <- curr_spearman$estimate
  volcano_data[i, 3] <- -log10(curr_spearman$p.value)
  
  curr_pearson <- cor.test((BBC_VARscores_pairs_mat_removed_totest %>% ungroup() %>%
                              filter(taxon_species == viruses[i]))$cord,
                           (BBC_VARscores_pairs_mat_removed_totest %>% ungroup() %>%
                              filter(taxon_species == viruses[i]))$child,
                           method = "pearson")
  
  volcano_data[i, 4] <- curr_pearson$estimate
  volcano_data[i, 5] <- -log10(curr_pearson$p.value)
  
  curr_wilcox <- wilcox.test((BBC_VARscores_pairs_mat_removed_totest %>% ungroup() %>%
                                filter(taxon_species == viruses[i]))$cord,
                             (BBC_VARscores_pairs_mat_removed_totest %>% ungroup() %>%
                                filter(taxon_species == viruses[i]))$child)
  
  volcano_data[i, 6] <- -log10(curr_wilcox$p.value)
  
  curr_ttest <- t.test((BBC_VARscores_pairs_mat_removed_totest %>% ungroup() %>%
                          filter(taxon_species == viruses[i]))$cord,
                       (BBC_VARscores_pairs_mat_removed_totest %>% ungroup() %>%
                          filter(taxon_species == viruses[i]))$child)
  
  volcano_data[i, 7] <- -log10(curr_ttest$p.value)
  
  volcano_data[i, 8] <- mean((BBC_VARscores_pairs_mat_removed_totest %>% ungroup() %>%
                                filter(taxon_species == viruses[i]))$cord,
                             na.rm = T) - 
    mean((BBC_VARscores_pairs_mat_removed_totest %>% ungroup() %>%
            filter(taxon_species == viruses[i]))$child,
         na.rm = T)
  
  volcano_data[i, 9] <- length(which((BBC_VARscores_pairs_mat_removed_totest %>% ungroup() %>%
                                        filter(taxon_species == viruses[i]))$cord > 1.936591))
  volcano_data[i, 10] <- length(which((BBC_VARscores_pairs_mat_removed_totest %>% ungroup() %>%
                                         filter(taxon_species == viruses[i]))$cord <= 1.936591))
  volcano_data[i, 11] <- length(which((BBC_VARscores_pairs_mat_removed_totest %>% ungroup() %>%
                                         filter(taxon_species == viruses[i]))$child > 1.936591))
  volcano_data[i, 12] <- length(which((BBC_VARscores_pairs_mat_removed_totest %>% ungroup() %>%
                                         filter(taxon_species == viruses[i]))$child <= 1.936591))
  
  curr_fishtest <- fisher.test(matrix(data = c(volcano_data[i, 11], volcano_data[i, 12],
                                               volcano_data[i, 9], volcano_data[i, 10]), 
                                      nrow = 2, ncol = 2))
  
  volcano_data[i, 13] <- curr_fishtest$estimate
  volcano_data[i, 14] <- -log10(curr_fishtest$p.value)
  
  
  print(i)
  
}

volcano_data <- volcano_data %>% filter(virus %nin%
                                          c("Bovine coronavirus", "Cercopithecine herpesvirus 1", "Cercopithecine herpesvirus 16",
                                            "Human respiratory syncytial virus", "NA")) # removing animal and unspecified

write.csv(volcano_data, "/BBC/20240424 - BBC cord child volcano.csv")

volcano_data <- read.csv("BBC/Data/20240424 - BBC cord child volcano.csv")

#fisher volcano
names(volcano_data) <- names(volcano_data) %>% str_replace_all("\\.", " ")

volcano_data <- volcano_data %>% mutate(`-log10(fisher fdr)` = -log10(p.adjust(10^-`fisher p`))) %>%
  mutate(`log2(odds ratio)` = log2(`odds ratio`))

curr_plot = ggplot() + geom_point(data = volcano_data, #%>% mutate(`odds ratio` = ifelse(`odds ratio` < 2^-15, 2^-15, `odds ratio`)),
                                  aes(x = `log2(odds ratio)`,
                                      y = `-log10(fisher fdr)`),
                                  alpha = .4) +
  theme_light() +
  theme(text = element_text(size = 8))

ggsave("BBC/Figures/20250324 - cord child fisher test volcano.jpeg",
       curr_plot,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300)

volcano_data <- volcano_data %>% mutate(virus = str_remove_all(virus, "Human "))

volcano_data_toplot <- volcano_data %>% filter(virus %nin% c("mastadenovirus C", "Norovirus MD145")) %>% #remove species duplicates
  mutate(virus = case_when(
    virus == "respiratory syncytial virus A" ~ "RSV A",
    virus == "respiratory syncytial virus B" ~ "RSV B",
    virus == "Epstein-Barr virus" ~ "EBV",
    virus == "herpesvirus 1" ~ "HSV 1",
    virus == "herpesvirus 2" ~ "HSV 2",
    virus == "herpesvirus 3" ~ "VZV",
    virus == "herpesvirus 6B" ~ "HHV 6B",
    virus == "adenovirus C serotype 2" ~ "Adeno C",
    virus == "adenovirus B serotype 16" ~ "Adeno B",
    virus == "Influenza B virus" ~ "Flu B",
    virus == "cytomegalovirus" ~ "CMV",
    virus == "parechovirus 1" ~ "pEV1",
    virus == "parechovirus 2" ~ "pEV2",
    virus == "Rhinovirus A" ~ "RV A",
    virus == "Rhinovirus B" ~ "RV B",
    virus == "Enterovirus A" ~ "EV A",
    virus == "Enterovirus B" ~ "EV B",
    virus == "Enterovirus C" ~ "EV C",
    virus == "Enterovirus D" ~ "EV D",
    virus == "Lordsdale virus" ~ "Lordsdale",
    T ~ virus
  ))

# top 11 differentially reactive viruses
f1 = ggplot() + geom_col(data = volcano_data_toplot %>% filter(`fisher p` > 75) %>%
                           mutate(`child` = `child pos` / (`child pos` + `child neg`),
                                  `cord` = `cord pos` / (`cord pos` + `cord neg`)) %>%
                           dplyr::select(virus, `child`, `cord`, `odds ratio`) %>%
                           pivot_longer(cols = c(`cord`, `child`), names_to = "Antibody Source",
                                        values_to = "Fraction Reactive") %>%
                           arrange(`odds ratio`),
                         aes(x = `Fraction Reactive`, y = virus, color = `Antibody Source`,
                             fill = `Antibody Source`), alpha = .5,
                         position = position_dodge2(padding = 0.2)) +
  scale_y_discrete(limits = rev((volcano_data_toplot %>% filter(`fisher p` > 75) %>% arrange(`odds ratio`))$virus)) +
  theme_light() +
  scale_fill_manual(values = c("black", "red")) +
  scale_color_manual(values = c("black", "red")) +
  theme(text = element_text(size = 8)) +
  theme(legend.position = "none", axis.title = element_blank())

# top 13 reactive viruses in child
f2 = ggplot() + geom_col(data = volcano_data_toplot %>% filter(`child pos` > 400) %>%
                           mutate(`child` = `child pos` / (`child pos` + `child neg`),
                                  `cord` = `cord pos` / (`cord pos` + `cord neg`)) %>%
                           dplyr::select(virus, `child`, `cord`, `odds ratio`) %>%
                           pivot_longer(cols = c(`cord`, `child`), names_to = "Antibody Source",
                                        values_to = "Fraction Reactive") %>%
                           arrange(`odds ratio`),
                         aes(x = `Fraction Reactive`, y = virus, color = `Antibody Source`,
                             fill = `Antibody Source`), alpha = .5,
                         position = position_dodge2(padding = 0.2)) +
  scale_y_discrete(limits = ((volcano_data_toplot %>% filter(`child pos` > 400) %>% arrange(`child pos`))$virus)) +
  theme_light() +
  scale_fill_manual(values = c("black", "red")) +
  scale_color_manual(values = c("black", "red")) +
  theme(text = element_text(size = 8)) +
  theme(legend.position = "none", axis.title = element_blank())


# top 15 reactive viruses in cord
f3 = ggplot() + geom_col(data = volcano_data_toplot %>% filter(`cord pos` > 650) %>%
                           mutate(`child` = `child pos` / (`child pos` + `child neg`),
                                  `cord` = `cord pos` / (`cord pos` + `cord neg`)) %>%
                           dplyr::select(virus, `child`, `cord`, `odds ratio`) %>%
                           pivot_longer(cols = c(`cord`, `child`), names_to = "Antibody Source",
                                        values_to = "Fraction Reactive") %>%
                           arrange(`odds ratio`),
                         aes(x = `Fraction Reactive`, y = virus, color = `Antibody Source`,
                             fill = `Antibody Source`), alpha = .5,
                         position = position_dodge2(padding = 0.2)) +
  scale_y_discrete(limits = ((volcano_data_toplot %>% filter(`cord pos` > 650) %>% arrange(`cord pos`))$virus)) +
  theme_light() +
  scale_fill_manual(values = c("black", "red")) +
  scale_color_manual(values = c("black", "red")) +
  theme(text = element_text(size = 8)) +
  theme(axis.title = element_blank())

g1 <- ggplotGrob(f1)
g2 <- ggplotGrob(f2)
g3 <- ggplotGrob(f3)

g <- cbind(g1, g2, g3)

grid.newpage()
grid.draw(g)

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - top differentially reactive cord child.jpeg",
       g,
       scale = 1,
       height = 2,
       width = 7.5,
       units = "in",
       dpi = 300)


##############
# are we detecting maternal antibodies in kid serum
##############

phiphfc_1 <- fread("BBC/Data/phipseq_0232_0241_VirscanLar_000_Hits_foldchange_annotated.tsv") %>%
  filter(taxon_species %in% BBC_VARscores_pairs_mat_removed_totest$taxon_species)
phiphfc_2 <- fread("BBC/Data/phipseq_0242_0254_VirscanLar_000_Hits_foldchange_annotated.tsv") %>%
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

sample_meta <- sample_meta %>% filter(lab_id %in% names(BBC_hfc))

IDs_torun <- sample_meta %>% group_by(ID) %>% summarise(n = n()) %>%
  filter(n == 2)

sample_meta <- sample_meta %>% filter(ID %in% IDs_torun$ID)

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
correlation_data <- fread("BBC/Data/20240315 - BBC cord child correlation.csv")

#pearson by age
corr_plot = ggplot(data = correlation_data,
                   aes(y = `Spearman rho`,
                       x = `Child Age`)) +
  geom_point(alpha = .5, size = .5) +
  theme_light() +
  geom_smooth() +
  theme(text = element_text(size = 8))

ggsave("BBC/Figures/20250324 - cord child spearman vs age.jpeg",
       corr_plot,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300) 