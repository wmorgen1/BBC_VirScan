
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

#######
# Data Prep
#######
#phipfc <- fread("C:/Users/wmorgen1/Documents/MARCC/BBC/phipseq_0232_0241_VirscanLar_000_FoldChange_annotated.tsv")
#phiphfc <- fread("C:/Users/wmorgen1/Documents/MARCC/BBC/phipseq_0232_0241_VirscanLar_000_Hits_foldchange_annotated.tsv")
phipfc <- fread("C:/Users/wmorgen1/Documents/MARCC/BBC/phipseq_0242_0254_VirscanLar_000_FoldChange_annotated.tsv")
phiphfc <- fread("C:/Users/wmorgen1/Documents/MARCC/BBC/phipseq_0242_0254_VirscanLar_000_Hits_foldchange_annotated.tsv")

entero_labels <- fread("C:/Users/wmorgen1/Documents/R/Larman_Lab/Helmsley/entero_rhino_cosa_species_label.csv")

phipfc <- phipfc %>% left_join(entero_labels) %>% 
  mutate(taxon_species=ifelse(taxon_species %in% entero_labels$taxon_species, ViralSpecies_label, taxon_species)) %>% 
  dplyr::select(-ViralSpecies_label) 

phiphfc <- phiphfc %>% left_join(entero_labels) %>% 
  mutate(taxon_species=ifelse(taxon_species %in% entero_labels$taxon_species, ViralSpecies_label, taxon_species)) %>% 
  dplyr::select(-ViralSpecies_label) 

phipfc <- phipfc %>% dplyr::select(
  "u_pep_id", "pep_id", "pos_start", "pos_end", "UniProt_acc", "pep_aa",       
  "taxon_genus", "taxon_species", "gene_symbol", "product",
  contains("beads"), contains("Beads"), contains("BEADS"),
  starts_with("252"), starts_with("253"), starts_with("254")
)
phiphfc <- phiphfc %>% dplyr::select(
  "u_pep_id", "pep_id", "pos_start", "pos_end", "UniProt_acc", "pep_aa",       
  "taxon_genus", "taxon_species", "gene_symbol", "product",
  contains("beads"), contains("Beads"), contains("BEADS"),
  starts_with("252"), starts_with("253"), starts_with("254")
)

################
# calculate Standard ARscores
################
# doubled the number of values of N for which random distributions are generated - 20240129 - reduced back to original
# instead of interpolating gamma parameters with a linear distribution, use a smoothed spline
# runs much slower

calc_scores <- function(norm_log, all_peptide_fcs, positives, positive_peps = NULL) {
  print("running ARscore algorithm")
  
  #representations <- c(
  #  30, 42, 60, 85, 120, 170, 240,
  #  340, 480, 679, 960, 1358, 1960,
  #  2715, 3920, 5431, 7840
  #  #, 2715, 3920, 5431, 7840 not needed for monkey virus, toxome
  #  #, 15680 only needed for arboscan
  #)
  representations <- c(
    30, 60, 120, 240,
    480, 960, 1960,
    3920, 7840
    #, 2715, 3920, 5431, 7840 not needed for monkey virus, toxome
    #, 15680 only needed for arboscan
  )  
  
  all_peptide_fcs <- all_peptide_fcs %>% mutate(fc = (value))
  
  # Option to change this to filtering by peptide
  all_peptide_fcs <- all_peptide_fcs %>% 
    filter(taxon_genus %nin% positives$taxon_genus)
  
  
  distributions <- data.frame(matrix(nrow = length(representations), ncol = 1000))
  row.names(distributions) <- representations
  
  #calculate random distributions
  set.seed(120)
  for(i in c(1:length(representations))){
    n = row.names(distributions)[i]
    for(j in c(1:1000)){
      distributions[i, j] = mean(sample(all_peptide_fcs$fc, as.numeric(n), replace = F))
      if(distributions[i, j] == 0) {distributions[i,j] = .001} 
    }
    distributions[i, 1000] = mean(all_peptide_fcs$fc)
    #print("calculating distributions " %+% i %+% " of " %+% length(representations))
  }
  
  #model random distributions as a gamma distribution
  dist_info <- distributions[,c(1, 2, 3)]
  names(dist_info) <- c("total_peps", "shape", "rate")
  for(i in c(1:length(representations))){
    
    dist_info[i, 1] <- row.names(dist_info)[i] %>% as.numeric()
    
    gammafit  <-  fitdist(distributions[i,] %>% t() %>% as.numeric(), "gamma")
    
    dist_info[i, 2] <- gammafit[["estimate"]][["shape"]]
    dist_info[i, 3] <- gammafit[["estimate"]][["rate"]]
  }
  
  dist_info <- dist_info %>% mutate(mean = shape / rate) %>%
    mutate(variance = shape / rate^2)
  
  # gamma distribution parameters vary with n
  #shape_lm <- lm(log10(shape) ~ I(log10(total_peps)), data = dist_info)
  shape_spline <- smooth.spline(x = log10(dist_info$total_peps), y = log10(dist_info$shape), spar = .5)
  #rate_lm <- lm(log10(rate) ~ I(log10(total_peps)), data = dist_info)
  rate_spline <- smooth.spline(x = log10(dist_info$total_peps), y = log10(dist_info$rate), spar = .5)
  #mean_lm <- lm(log10(mean) ~ I(log10(total_peps)), data = dist_info)
  mean_spline <- smooth.spline(x = log10(dist_info$total_peps), y = log10(dist_info$mean), spar = .5)
  
  # add distribution info to norm_log
  norm_log <- norm_log %>%
    mutate(shape = 10^(predict(shape_spline, log10(total_peps))$y),
           rate = 10^(predict(rate_spline, log10(total_peps))$y)) %>%
    mutate(mean = shape / rate) %>%
    mutate(variance = shape / rate^2)
  
  # calculate scores and p values
  norm_log <- norm_log %>%
    mutate(vir_score = zscoreGamma(score_norm, shape = shape / sqrt(total_peps), rate = rate / sqrt(total_peps)),
           virus_fc = score_norm / mean) %>%
    mutate(vir_score = ifelse(vir_score < -10, -10, vir_score)) %>%
    mutate(p_val = -log10((1-pnorm(abs(zscoreGamma(score_norm, shape = shape, rate = rate)))))) %>%
    mutate(p_val = ifelse(p_val>15, 15, p_val))
  
  return(norm_log)
}

iterative_scores <- function(norm_log_1, all_peptide_fcs_1, max_iterations = 10,
                             p_cutoff = -log10(.0001), score_cutoff = 2) {
  iterations = 0
  positives_1 = norm_log_1 %>% filter(F)
  positives_output = list()
  library(limma)
  
  while(iterations < max_iterations) {
    print("iteration " %+% (iterations+1))
    scores <- calc_scores(norm_log = norm_log_1, all_peptide_fcs = all_peptide_fcs_1, positives = positives_1)
    
    #update variables
    positives_2 <- scores %>% filter(p_val == 15) %>% filter(vir_score > 0) #12/02/23 added second internal criteria
    positives_1 <- scores %>% filter(p_val > p_cutoff) %>% filter(vir_score > score_cutoff) %>% full_join(positives_2, 
                         join_by(taxon_species, taxon_genus, sample_id, total_peps, score, score_norm, shape, 
                         rate, mean, variance, vir_score, virus_fc, p_val))
    positives_output[[iterations+1]] <- positives_1 
    
    #add break if no new positives
    if(iterations > 0) {
      if((positives_1$taxon_species %>% length()) <= (positives_output[[iterations]]$taxon_species %>% length())) {
        iterations <- iterations + 10
      }
    }
    
    iterations <- iterations + 1
  }
  
  outputs <- list()
  outputs[[1]] <- scores
  outputs[[2]] <- positives_output
  return(outputs) 
}
# added sublibrary argument
ARscore_algorithm <- function(hfc = NULL, fc, set_max_iterations = 10,
                              set_p_cutoff = -log10(.0001), set_score_cutoff = 2,
                              required_number_of_peptides = 50, sublib = NA) {
  
  ## peptides that have hits in beads
  if(!is.null(hfc)){
    bad_beads <- hfc %>% dplyr::select(pep_aa,contains('BEADS'),contains('Beads'),contains('beads')) %>% pivot_longer(-pep_aa) %>%
      filter(value > 1) %>% dplyr::select(pep_aa) %>% pull
    
    ### fold change ###
    d1 <- fc
    
    ## make longform dataframe out of fc, no longer floored fc at 1
    d2 <- d1 %>% filter(pep_aa %nin% bad_beads) %>% 
      dplyr::select(-c(contains('BEADS'),contains('Beads'),contains('beads'))) %>%
      #common pulldowns. Has to be edited in some screens
      pivot_longer(cols = c(contains('20A20G'), contains('20S'),
                            contains('ugIg'), contains('DPI'),
                            contains('IgA'), contains('-')),names_to = 'sample_id')
    
  } else {
    
    ### fold change ###
    d1 <- fc
    
    ## make longform dataframe out of fc, no longer floored fc at 1
    d2 <- d1 %>% 
      dplyr::select(-c(contains('BEADS'),contains('Beads'),contains('beads'))) %>%
      #common pulldowns. Has to be edited in some screens
      pivot_longer(cols = c(contains('20A20G'), contains('20S'),
                            contains('ugIg'), contains('DPI'),
                            contains('IgA'), contains('-')),names_to = 'sample_id')
    
  }
 
  ids <- d2 %>% dplyr::select(sample_id) %>% unique()
  
  # need to group by taxon_species, sample_id and count number of reactivities. Then add back sample annotations
  if (is.na(sublib)) {
    d4 <- d2 %>% group_by(taxon_species, taxon_genus, sample_id) %>% summarise(score = sum((value))) %>% left_join(ids)
    d4_2 <- d2 %>% group_by(taxon_species, taxon_genus, sample_id) %>% summarise(total_peps = n()) %>%
      left_join(ids) %>% left_join(d4) %>% mutate(score = ifelse(is.na(score), 0, score))
    d4_2 <- d4_2 %>% mutate(score_norm = score / total_peps) %>% filter(total_peps >= required_number_of_peptides)
  } else {
    ### adding grouping here the new phageome library has several sublibraries: c("pepsyn","dolphin","dolphinepitopes")
    ## one or more of these can be included. They are added as "sublib" in annotations
    d4 <- d2 %>% filter(sublibrary %in% sublib) %>% group_by(taxon_species, taxon_genus, sample_id) %>% summarise(score = sum((value))) %>% left_join(ids)
    d4_2 <- d2 %>% filter(sublibrary %in% sublib) %>% group_by(taxon_species, taxon_genus, sample_id) %>% summarise(total_peps = n()) %>%
      left_join(ids) %>% left_join(d4) %>% mutate(score = ifelse(is.na(score), 0, score))
    d4_2 <- d4_2 %>% mutate(score_norm = score / total_peps) %>% filter(total_peps >= required_number_of_peptides)
    
  }
  
  library(tidyverse)
  library(fitdistrplus)
  library(limma)
  
  algorithm_output <- list()
  
  for(R in c(1:length(ids$sample_id))) {
    print(R)
    algorithm_output[[R]] <- iterative_scores(norm_log_1 = d4_2 %>% filter(sample_id == ids$sample_id[R]), 
                                              all_peptide_fcs_1 = d2 %>% filter(sample_id == ids$sample_id[R]),
                                              max_iterations = set_max_iterations,
                                              p_cutoff = set_p_cutoff, 
                                              score_cutoff = set_score_cutoff)
    
  }
  
  #algorithm output has the final output as well as a structure that tracks positives 
  scores <- algorithm_output[[1]][[1]]
  if(length(algorithm_output)>1){
    for(i in c(2:length(algorithm_output))){
      scores <- rbind(scores, algorithm_output[[i]][[1]])
    }
  }
  
  return(scores)
}

phip_VARscores <- ARscore_algorithm(fc = phipfc, hfc = phiphfc,
                                    required_number_of_peptides = 50)

write_rds(phip_VARscores, "C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 11.rds")

#############
# make wide (if many samples were run)
#############

VARscores_1 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 1.rds")
VARscores_2 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 2.rds")
VARscores_3 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 3.rds")
VARscores_4 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 4.rds")
VARscores_5 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 5.rds")
VARscores_6 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 6.rds")
VARscores_7 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 7.rds")
VARscores_8 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 8.rds")
VARscores_9 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 9.rds")
VARscores_10 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 10.rds")
VARscores_11 <- read_rds("C:/Users/wmorgen1/Documents/MARCC/BBC/20240206 - VARscores - 11.rds")

VARscores_wide <- VARscores_1 %>% full_join(VARscores_2) %>% 
  full_join(VARscores_3) %>% full_join(VARscores_4) %>%
  full_join(VARscores_5) %>% full_join(VARscores_6) %>%
  full_join(VARscores_7) %>% full_join(VARscores_8) %>%
  full_join(VARscores_9) %>% full_join(VARscores_10) %>%
  full_join(VARscores_11) %>%
  dplyr::select(taxon_species, taxon_genus, total_peps, sample_id, vir_score) %>%
  pivot_wider(names_from = sample_id, values_from = vir_score)

write.csv(VARscores_wide, "C:/Users/wmorgen1/Documents/MARCC/BBC/20240214 - VRC VARscores.csv")

#####
ebv_cmv_HSV_rsv_HIV <- VARscores_1 %>% full_join(VARscores_2) %>% 
  full_join(VARscores_3) %>% full_join(VARscores_4) %>%
  full_join(VARscores_5) %>% full_join(VARscores_6) %>%
  full_join(VARscores_7) %>% full_join(VARscores_8) %>%
  full_join(VARscores_9) %>% full_join(VARscores_10) %>%
  filter(taxon_species %in% c("Human cytomegalovirus", "Epstein-Barr virus",
                              "Human immunodeficiency virus type 1 group M subtype C", 
                              "Human respiratory syncytial virus A",
                              "Hepatitis C virus genotype 1b", "Zaire ebolavirus",
                              "Enterovirus B", 
                              "Human herpesvirus 1",
                              "Reston ebolavirus")) %>%
  mutate(taxon_species = as.factor(taxon_species)) %>%
  mutate(taxon_species = fct_relevel(taxon_species, "Human respiratory syncytial virus A", 
                                     "Enterovirus B", "Epstein-Barr virus",
                                     "Human cytomegalovirus", "Human herpesvirus 1",
                                     "Human immunodeficiency virus type 1 group M subtype C",
                                     "Hepatitis C virus genotype 1b", "Zaire ebolavirus",
                                     "Reston ebolavirus")) %>%
  mutate(vir_score = ifelse(vir_score>30, 30, vir_score))

ggplot() + 
  geom_point(data = ebv_cmv_HSV_rsv_HIV, aes(x = taxon_species, y = vir_score),
             position = position_quasirandom(),
             size = 1,
             alpha = .3) +
  scale_x_discrete(labels=c("Human respiratory syncytial virus A" = "RSV",
                            "Enterovirus B" = "Entero B",
                            "Epstein-Barr virus" = "EBV",
                            "Human cytomegalovirus" = "CMV",
                            "Human herpesvirus 1" = "HSV1",
                            "Human immunodeficiency virus type 1 group M subtype C" = "HIV",
                            "Hepatitis C virus genotype 1b" = "HCV",
                            "Zaire ebolavirus" = "Ebola",
                            "Reston ebolavirus" = "Ebola 2")) +
  theme_classic() +
  geom_hline(yintercept = 2)

# HSV1 HSV2 comparison

HSV1_HSV2 <- VARscores_1 %>% full_join(VARscores_2) %>% 
  full_join(VARscores_3) %>% full_join(VARscores_4) %>%
  full_join(VARscores_5) %>% full_join(VARscores_6) %>%
  full_join(VARscores_7) %>% full_join(VARscores_8) %>%
  full_join(VARscores_9) %>% full_join(VARscores_10) %>% ungroup() %>%
  filter(taxon_genus == "Simplexvirus") %>%
  pivot_wider(id_cols = sample_id, names_from = taxon_species,
              values_from = vir_score, values_fn = max)

ggplot() + 
  geom_point(data = HSV1_HSV2,
             aes(x = `Human herpesvirus 1`, 
                 y = `Human herpesvirus 2`),
             size = 1,
             alpha = .2) +
  xlab("HSV1") +
  ylab("HSV2") +
  geom_hline(yintercept = 2) +
  geom_vline(xintercept = 2) +
  theme_light()

#####

HSV1_HSV2 <- VARscores_1 %>% full_join(VARscores_2) %>% 
  full_join(VARscores_3) %>% full_join(VARscores_4) %>%
  full_join(VARscores_5) %>% full_join(VARscores_6) %>%
  full_join(VARscores_7) %>% full_join(VARscores_8) %>%
  full_join(VARscores_9) %>% full_join(VARscores_10) %>% ungroup() %>%
  filter(taxon_species %in% c("Reston ebolavirus",
                              "Human respiratory syncytial virus A")) %>%
  pivot_wider(id_cols = sample_id, names_from = taxon_species,
              values_from = vir_score, values_fn = max)

ggplot() + 
  geom_point(data = HSV1_HSV2,
             aes(x = `Reston ebolavirus`, 
                 y = `Human respiratory syncytial virus A`),
             size = 1,
             alpha = .2) +
  xlab("Ebola") +
  ylab("RSV A") +
  geom_hline(yintercept = 2) +
  geom_vline(xintercept = 2) +
  theme_light()

HSV1_HSV2 <- VARscores_1 %>% full_join(VARscores_2) %>% 
  full_join(VARscores_3) %>% full_join(VARscores_4) %>%
  full_join(VARscores_5) %>% full_join(VARscores_6) %>%
  full_join(VARscores_7) %>% full_join(VARscores_8) %>%
  full_join(VARscores_9) %>% full_join(VARscores_10) %>% ungroup() %>%
  filter(taxon_species %in% c("Reston ebolavirus",
                              "Human cytomegalovirus")) %>%
  pivot_wider(id_cols = sample_id, names_from = taxon_species,
              values_from = vir_score, values_fn = max)

ggplot() + 
  geom_point(data = HSV1_HSV2,
             aes(x = `Reston ebolavirus`, 
                 y = `Human cytomegalovirus`),
             size = 1,
             alpha = .2) +
  xlab("Ebola") +
  ylab("CMV") +
  geom_hline(yintercept = 2) +
  geom_vline(xintercept = 2) +
  theme_light()
#####