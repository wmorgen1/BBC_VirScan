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

# excluding unspecified
hsv_VAR <- BBC_VARscores %>%
  filter(taxon_genus %in% c("Orthopneumovirus", "Metapneumovirus")) %>%
  filter(taxon_species != "Human respiratory syncytial virus") %>%
  mutate(taxon_species = as.factor(taxon_species)) %>%
  arrange(taxon_species)


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
        legend.position = "none") + 
  geom_hline(yintercept = 1.936591)+ 
  geom_vline(xintercept = 1.936591)+
  geom_abline(intercept = 0, slope = 2)+
  geom_abline(intercept = 0, slope = 0.5)

######
# peptide epitopes
######
child_meta <- sample_meta %>% filter(Plasma_type =="child")

phip1 <- fread("C:/Users/wmorgen1/Documents/MARCC/BBC/phipseq_0242_0254_VirscanLar_000_Hits_foldchange_annotated.tsv") %>%
  filter(taxon_species %in% c("Human respiratory syncytial virus A", "Human respiratory syncytial virus B"))
phip2 <- fread("C:/Users/wmorgen1/Documents/MARCC/BBC/phipseq_0232_0241_VirscanLar_000_Hits_foldchange_annotated.tsv") %>%
  filter(taxon_species %in% c("Human respiratory syncytial virus A", "Human respiratory syncytial virus B"))

phip <- phip1 %>% full_join(phip2)

names(phip) <- names(phip) %>% str_replace_all("-", ".")

child_phip <- phip %>%
  dplyr::select(c("pep_id", "pos_start", "pos_end",
                  "UniProt_acc", "pep_aa", "taxon_genus",
                  "taxon_species", "gene_symbol", "product", child_meta$lab_id %>% str_remove("^X")))

child_phip_long <- child_phip %>% pivot_longer(cols = contains("Xi"), names_to = "lab_id") %>%
  mutate(lab_id = "X" %+% lab_id) 

hits <- child_phip_long %>% group_by(pep_id) %>% 
  summarise(hits = sum(value>1), n = n()) %>% filter(hits > 9)

child_phip_long <- child_phip_long %>% filter(pep_id %in% hits$pep_id)

RSV_VAR_forcors <- hsv_VAR %>% dplyr::select(taxon_species, lab_id, vir_score) %>%
  pivot_wider(names_from = taxon_species, values_from = vir_score)

child_phip_long <- child_phip_long %>% left_join(RSV_VAR_forcors)

pep_VAR_correlations <- hits %>% mutate(`RSV A Pearson R` = 0, `RSV A p value` = 0,
                                        `RSV B Pearson R` = 0, `RSV B p value` = 0)

for(i in c(1:length(hits$pep_id))){
  curr_pep <- pep_VAR_correlations$pep_id[i]
  
  curr_cor <- cor.test((child_phip_long %>% filter(pep_id == curr_pep))$value,
                       (child_phip_long %>% filter(pep_id == curr_pep))$`Human respiratory syncytial virus A`)
  
  pep_VAR_correlations$`RSV A Pearson R`[i] = curr_cor$estimate
  pep_VAR_correlations$`RSV A p value`[i] = -log10(curr_cor$p.value)
  
  curr_cor <- cor.test((child_phip_long %>% filter(pep_id == curr_pep))$value,
                       (child_phip_long %>% filter(pep_id == curr_pep))$`Human respiratory syncytial virus B`)
  
  pep_VAR_correlations$`RSV B Pearson R`[i] = curr_cor$estimate
  pep_VAR_correlations$`RSV B p value`[i] = -log10(curr_cor$p.value)
  
}

write_rds(pep_VAR_correlations, "C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/20240503 - RSV pep VAR cors.rds")

pep_meta <- phip %>% dplyr::select("pep_id"               ,                                               
                                   "pos_start"            ,                                              
                                   "pos_end"              ,                                             
                                   "UniProt_acc"          ,                                            
                                   "pep_aa"               ,                                           
                                   "taxon_genus"          ,                                          
                                   "taxon_species"        ,                                         
                                   "gene_symbol"          ,                                        
                                   "product")

pep_VAR_correlations <- pep_VAR_correlations %>%
  left_join(pep_meta)

ggplot() + geom_point(data = pep_VAR_correlations %>% 
                        filter(gene_symbol != "NA") %>%
                        mutate(Protein = gene_symbol),
                      aes(x = `RSV A Pearson R`,
                      y = `RSV B Pearson R`,
                      color = Protein),
                      size = .6,
                      alpha = .5) +
  theme_light() + 
  xlab("Correlation Coefficient\nBetween Peptide Reactivity\nand RSV A VARscore") +
  ylab("Correlation Coefficient\nBetween Peptide Reactivity\nand RSV B VARscore") +
  scale_color_manual(values = c(
    "#1E88E5",
    "#D81B60",
    "#004D40",
    "#FFC107",
    "grey70"
  )) + 
  theme(text = element_text(size = 8))

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240503 - RSV pep by VAR cors - color by protein.jpeg",
       scale = 1,
       height = 2.5,
       width = 3.5,
       units = "in",
       dpi = 300)


ggplot() + geom_point(data = pep_VAR_correlations %>% filter(gene_symbol != "NA") %>%
                        mutate(Type = ifelse(taxon_species == "Human respiratory syncytial virus A", "A", "B")),
                      aes(x = `RSV A Pearson R`,
                          y = `RSV B Pearson R`,
                          color = Type),
                      size = .6,
                      alpha = .5) +
  theme_light() + 
  xlab("Correlation Coefficient\nBetween Peptide Reactivity\nand RSV A VARscore") +
  ylab("Correlation Coefficient\nBetween Peptide Reactivity\nand RSV B VARscore") +
  scale_color_manual(values = c("#005AB5", "#DC3220")) + 
  theme(text = element_text(size = 8))

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240503 - RSV pep by VAR cors - color by type.jpeg",
       scale = 1,
       height = 2.5,
       width = 3,
       units = "in",
       dpi = 300)


# GENOME PLOT FOR G PROTEIN
library(Biostrings)

pep_VAR_correlations_genomeplot <- pep_VAR_correlations %>% filter(gene_symbol == "G")

ref_g_seq <- readAAStringSet("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/protein_sequences/RSV S2 ts1C Glycoprotein G.fasta")
#need to align peptides to the reference protein

rsv_g_alignments <- pairwiseAlignment(pep_VAR_correlations_genomeplot$pep_aa, ref_g_seq, type = "global-local")

pep_VAR_correlations_genomeplot <- pep_VAR_correlations_genomeplot %>%
  mutate(start = rsv_g_alignments@subject@range@start) %>%
  mutate(end = rsv_g_alignments@subject@range@start + rsv_g_alignments@subject@range@width)

# annotations of the protein from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4211642/#:~:text=The%20RSV%20A2%20strain%20G,a%20pair%20of%20disulfide%20bonds.
# 1-36 
# 36-67 Transmembrane
# 67-164 Mucin-like region I
# 164-186 Central conserved domain with cysteine noose
# 186-224 Heparin Binding Domain
# 224-298 Mucin-like region II

RSV_G_labels <- data.frame(matrix(nrow = 6))
RSV_G_labels$this_row <- c(1,2,3,4,5,6)
RSV_G_labels$subsegment <- c("", "TM", "MLR I", "CCD", "HBD", "MLR II")
RSV_G_labels$start <- c(1,36,67,164,186,224)
RSV_G_labels$stop <- c(36,67,164,186,224, 298)

RSV_G <- ggplot() +
  geom_hline(yintercept = 0, color = "black") +
  geom_rect(data = pep_VAR_correlations_genomeplot %>%
              mutate(Type = ifelse(taxon_species == "Human respiratory syncytial virus A", "A", "B")),
            aes(xmin = start,
                xmax = end,
                ymin = `RSV A Pearson R` - 0.03,
                ymax = `RSV A Pearson R`,
                fill = Type),
            alpha = .6,
            color = "black") +
  scale_fill_manual(values = c("#005AB5", "#DC3220")) +
  ylab("Correlation Coefficient\nBetween Peptide Reactivity\nand RSV A VARscore") +
  scale_y_continuous(limits = c(-.25, 0.95)) + 
  scale_x_continuous(breaks = c(36, 67, 164, 186, 224),
                     limits = c(0 , 299)) + 
  xlab("Aligned RSV G position") + theme_light() + 
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        text = element_text(size=8), 
        axis.text.x = element_blank())  + 
  #X-axis label here with genome
  geom_rect(data = RSV_G_labels,
            aes(xmin=start, xmax = stop,
                ymin=if_else(this_row %%2 == 0,-.14, -.16), 
                ymax=if_else(this_row %%2 == 0,-.16, -.18)),
            color = "grey50",alpha = 0.7) +
  geom_text(data = RSV_G_labels, aes(x = (start+stop)/2,label = subsegment,y=-.21,angle = 0), size = 2)

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240503 - RSV G genome plot pep by VAR cors - color by type.jpeg",
       RSV_G,
       scale = 1,
       height = 2,
       width = 6,
       units = "in",
       dpi = 300)

RSV_G <- ggplot() +
  geom_hline(yintercept = 0, color = "black") +
  geom_rect(data = pep_VAR_correlations_genomeplot %>%
              mutate(Type = ifelse(taxon_species == "Human respiratory syncytial virus A", "A", "B")),
            aes(xmin = start,
                xmax = end,
                ymin = `RSV B Pearson R` - 0.03,
                ymax = `RSV B Pearson R`,
                fill = Type),
            alpha = .6,
            color = "black") +
  scale_fill_manual(values = c("#005AB5", "#DC3220")) +
  ylab("Correlation Coefficient\nBetween Peptide Reactivity\nand RSV B VARscore") +
  scale_y_continuous(limits = c(-.25, 0.95)) + 
  scale_x_continuous(breaks = c(36, 67, 164, 186, 224),
                     limits = c(0 , 299)) + 
  xlab("Aligned RSV G position") + theme_light() + 
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        text = element_text(size=8), 
        axis.text.x = element_blank())  + 
  #X-axis label here with genome
  geom_rect(data = RSV_G_labels,
            aes(xmin=start, xmax = stop,
                ymin=if_else(this_row %%2 == 0,-.14, -.16), 
                ymax=if_else(this_row %%2 == 0,-.16, -.18)),
            color = "grey50",alpha = 0.7) +
  geom_text(data = RSV_G_labels, aes(x = (start+stop)/2,label = subsegment,y=-.21,angle = 0), size = 2)


ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240503 - RSV G genome plot pep by B VAR cors - color by type.jpeg",
       RSV_G,
       scale = 1,
       height = 2,
       width = 6,
       units = "in",
       dpi = 300)