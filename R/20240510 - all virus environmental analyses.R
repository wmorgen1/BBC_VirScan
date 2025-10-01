setwd("C:/Users/wmorgen1/Documents/Larman_Lab/BBC/")
library(tidyverse)
library(RColorBrewer)
library(gtable)
library(grid)
library(ggplot2)
library(ggbeeswarm)
library(data.table)
library(ARscore)

########
# Data Prep
#######
library(readxl)
sample_meta <- read_excel("sample_meta.xls")

viruses_totest <- read.csv("Data/20240424 - BBC cord child volcano.csv")

sample_meta <- sample_meta %>% mutate(age = difftime(time1 = Plasma_Date, time2 = DOD, units = "days")) %>%
  mutate(years = (age %>% as.numeric())/365.25) %>% mutate(Plasma_type = ifelse(Plasma_type=="Postnatal","child",Plasma_type))

breast_fed_meta <- read_excel("breastfed_update ID with codebook.xlsx")
parity_meta <- read_excel("parity.xlsx")
care_meta <- read_csv("child care.csv") %>% filter((ID != 1574 | I_childcare == 1)) %>%
  filter((ID != 2097 | I_cc1begyr == 0)) %>%
  filter((ID != 3836 | I_cc1begyr == 3)) %>%
  dplyr::select("ID", "I_childcare", "I_cc1begyr", "I_cc1begmo", "I_cc1begday") %>%
  mutate(I_cc1begyr = case_when(
    !is.na(I_cc1begyr) ~ I_cc1begyr,
    !is.na(I_cc1begmo) ~ 0,
    !is.na(I_cc1begday) ~ 0,
    is.na(I_cc1begyr) ~ I_cc1begyr
  )) %>%
  mutate(I_cc1begmo = case_when(
    !is.na(I_cc1begmo) ~ I_cc1begmo,
    !is.na(I_cc1begyr) ~ 0,
    !is.na(I_cc1begday) ~ 0,
    is.na(I_cc1begmo) ~ I_cc1begmo
  )) %>%
  mutate(I_cc1begday = case_when(
    !is.na(I_cc1begday) ~ I_cc1begday,
    !is.na(I_cc1begmo) ~ 0,
    !is.na(I_cc1begyr) ~ 0,
    is.na(I_cc1begday) ~ I_cc1begday
  )) %>% mutate(begin_of_cc = I_cc1begyr + I_cc1begmo / 12 + I_cc1begday / 365.25) %>% unique()
additional_mom_meta <- read_excel("pheno_update2.xlsx")

sample_meta <- sample_meta %>% left_join(breast_fed_meta) %>%
  left_join(parity_meta) %>% left_join(care_meta) %>%
  left_join(additional_mom_meta)

BBC_VARscores <- read_rds("Data/20240424 - BBC VARscores - maternal abs removed.rds")

BBC_VARscores <- BBC_VARscores %>% dplyr::rename(vir_score = value) %>%
  left_join(sample_meta, multiple = "all")

hsv_VAR <- BBC_VARscores %>%
  filter(taxon_species %in% viruses_totest$virus) %>%
  mutate(taxon_species = as.factor(taxon_species)) %>%
  mutate(race_update = as.factor(race_update)) %>%
  mutate(race_update = fct_relevel(race_update, c("2", "1", "3", "4", "5", "6", "7", "10")))%>%
  mutate(birthyear = as.factor(year(DOD))) %>%
  mutate(birthyear = fct_relevel(birthyear, c("2007", "2000", "2001", "2002",
                                              "2003", "2004", "2005", "2006",
                                              "2008", "2009", "2010", "2011")))%>%
  mutate(birthmonth = as.factor(month(DOD))) %>%
  mutate(birthmonth = fct_relevel(birthmonth, c("9", "1", "2", "3",
                                              "4", "5", "6", "7",
                                              "8", "10", "11", "12"))) %>%
  mutate(race_update = ifelse(race_update == "5", "1",race_update))

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
  mutate(`Mom Entero positive` = (`Enterovirus A` >= 1.936591) | (`Enterovirus B` >= 1.936591) |
           (`Enterovirus C` >= 1.936591) | (`Enterovirus D` >= 1.936591) |
           (`Rhinovirus B` >= 1.936591) | (`Rhinovirus A` >= 1.936591)) %>%
  mutate(`Mom EBV positive` = (`Epstein-Barr virus` >= 1.936591)) %>%
  mutate(`Mom Adeno positive` = (`Human adenovirus 14`>= 1.936591 |                
                     `Human adenovirus 55`>= 1.936591 |                 `Human adenovirus A serotype 12`>= 1.936591 |     
                     `Human adenovirus A serotype 18`>= 1.936591 |      `Human adenovirus A serotype 31`>= 1.936591 |     
                     `Human adenovirus B serotype 16`>= 1.936591 |      `Human adenovirus B serotype 3`>= 1.936591 |      
                     `Human adenovirus B serotype 7`>= 1.936591 |       `Human adenovirus C serotype 1`>= 1.936591 |      
                     `Human adenovirus C serotype 2`>= 1.936591 |       `Human adenovirus C serotype 5`>= 1.936591 |      
                     `Human adenovirus D serotype 15`>= 1.936591 |      `Human adenovirus D serotype 17`>= 1.936591 |     
                     `Human adenovirus D37`>= 1.936591 |                `Human adenovirus E serotype 4`>= 1.936591 |      
                     `Human adenovirus F serotype 40`>= 1.936591 |      `Human adenovirus F serotype 41`>= 1.936591)) %>%
  mutate(`Mom Astro positive` = (`Human astrovirus-1`>= 1.936591 | `Human astrovirus-6`>= 1.936591 | `Human astrovirus-8`>= 1.936591)) %>%
  mutate(`Mom Parvo positive` = (`Human parvovirus B19` >= 1.936591)) %>%
  mutate(`Mom HSV6 positive` = (`Human herpesvirus 6B` >= 1.936591) | (`Human herpesvirus 6A` >= 1.936591)) %>%
  mutate(`Mom HSV7 positive` = (`Human herpesvirus 7` >= 1.936591)) %>%
  mutate(`Mom RSV positive` = (`Human respiratory syncytial virus A` >= 1.936591) | (`Human respiratory syncytial virus B` >= 1.936591)) %>%
  mutate(`Mom Flu positive` = (`Influenza A virus`>= 1.936591 | `Influenza B virus`>= 1.936591 | `Influenza C virus`>= 1.936591)) %>%
  mutate(`Mom Noro positive` = (`Lordsdale virus`>= 1.936591 | `Norovirus MD145`>= 1.936591 | `Norwalk virus`>= 1.936591 | `Norwalk-like virus`>= 1.936591)) %>%
  mutate(`Mom VZV positive` = (`Varicella-zoster virus` >= 1.936591) | (`Human herpesvirus 3` >= 1.936591)) %>%
  dplyr::select(`ID`, contains("Mom"), `breast_fed1`) %>%
  mutate(`Mom CMV positive` = as.factor(`Mom CMV positive`),
         `Mom HSV1 positive` = as.factor(`Mom HSV1 positive`),
         `Mom HSV2 positive` = as.factor(`Mom HSV2 positive`),
         `Mom Entero positive` = as.factor(`Mom Entero positive`),
         `Mom EBV positive` = as.factor(`Mom EBV positive`),
         `Mom Adeno positive` = as.factor(`Mom Adeno positive`),
         `Mom Astro positive` = as.factor(`Mom Astro positive`),
         `Mom Noro positive` = as.factor(`Mom Noro positive`),
         `Mom Parvo positive` = as.factor(`Mom Parvo positive`),
         `Mom HSV6 positive` = as.factor(`Mom HSV6 positive`),
         `Mom HSV7 positive` = as.factor(`Mom HSV7 positive`),
         `Mom RSV positive` = as.factor(`Mom RSV positive`),
         `Mom Flu positive` = as.factor(`Mom Flu positive`),
         `Mom VZV positive` = as.factor(`Mom VZV positive`)
  )

#####
# Regression and Survival Analyses
#####

viruses_totest <- viruses_totest %>% filter(virus != "Human adenovirus 14") %>% 
  filter(virus != "Human adenovirus 55") %>% 
  filter(virus != "Human adenovirus A serotype 18") %>% 
  filter(virus != "Human adenovirus A serotype 31") %>% 
  filter(virus != "Human adenovirus B serotype 3") %>% 
  filter(virus != "Human adenovirus B serotype 7") %>% 
  filter(virus != "Human mastadenovirus C") %>% 
  filter(virus != "Human adenovirus C serotype 5") %>% 
  filter(virus != "Human adenovirus D serotype 15") %>% 
  filter(virus != "Human adenovirus D serotype 17") %>% 
  filter(virus != "Human adenovirus F serotype 41") %>%
  filter(virus != "Human mastadenovirus E") %>%
  filter(virus != "Norwalk-like virus") %>%
  filter(virus != "Human astrovirus-6") %>%
  filter(virus != "Human Norovirus MD145")

survival_outputs <- data.frame(matrix(nrow = 49, ncol = 91)) 
names(survival_outputs) <- c("Virus", 
                             "Childcare Estimate", "Childcare P",
                             "Breast Estimate", "Breast P",
                             "Maternal Age Estimate", "Maternal Age P",
                             "Preterm Estimate", "Preterm P",
                             "Parity Estimate", "Parity P",
                             "Delivery Estimate", "Delivery P",
                             "Sex Estimate", "Sex P",
                             "Race 1 Estimate", "Race 1 P",
                             "Race 3 Estimate", "Race 3 P",
                             "Race 4 Estimate", "Race 4 P",
                             "Race 6 Estimate", "Race 6 P",
                             "Race 7 Estimate", "Race 7 P",
                             "Race 10 Estimate", "Race 10 P",
                             "Birth 2001 Estimate", "Birth 2001 P",
                             "Birth 2002 Estimate", "Birth 2002 P",
                             "Birth 2003 Estimate", "Birth 2003 P",
                             "Birth 2004 Estimate", "Birth 2004 P",
                             "Birth 2005 Estimate", "Birth 2005 P",
                             "Birth 2006 Estimate", "Birth 2006 P",
                             "Birth 2008 Estimate", "Birth 2008 P",
                             "Birth 2009 Estimate", "Birth 2009 P",
                             "Birth 2010 Estimate", "Birth 2010 P",
                             "Birth 2011 Estimate", "Birth 2011 P",
                             "Birth Jan Estimate", "Birth Jan P",
                             "Birth Feb Estimate", "Birth Feb P",
                             "Birth Mar Estimate", "Birth Mar P",
                             "Birth Apr Estimate", "Birth Apr P",
                             "Birth May Estimate", "Birth May P",
                             "Birth June Estimate", "Birth June P",
                             "Birth July Estimate", "Birth July P",
                             "Birth Aug Estimate", "Birth Aug P",
                             "Birth Oct Estimate", "Birth Oct P",
                             "Birth Nov Estimate", "Birth Nov P",
                             "Birth Dec Estimate", "Birth Dec P",
                             "Mom CMV Estimate", "Mom CMV P",
                             "Mom HSV1 Estimate", "Mom HSV1 P",
                             "Mom HSV2 Estimate", "Mom HSV2 P",
                             "Mom EBV Estimate", "Mom EBV P",
                             "Mom Parvo Estimate", "Mom Parvo P",
                             "Mom HSV6 Estimate", "Mom HSV6 P",
                             "Mom HSV7 Estimate", "Mom HSV7 P",
                             "Mom RSV Estimate", "Mom RSV P",
                             "Mom VZV Estimate", "Mom VZV P",
                             "Age Estimate", "Age P"
)

for(i in c(1:49)){
  curr_virus <- viruses_totest$virus[i]
  
  survival_outputs$Virus[i] <- curr_virus
  
  curr_model_data <- hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
    left_join(Maternal_pos) %>% filter(taxon_species == curr_virus) %>%
    mutate(`kid pos` = as.factor(vir_score > 1.936591)) %>%
    mutate(cc = case_when(
      (years < begin_of_cc) | (I_childcare == 0) ~ "No Childcare",
      years > begin_of_cc ~ "Childcare",
      T ~ NA
    )) %>%
    mutate(breast_fed = as.factor(breast_fed1 %in% c(1, 2))) 
  
  curr_model <-  glm(`kid pos` ~ `cc` + `breast_fed` + mage + preterm1 + 
                       parity2 + DELTYPE + sex + race_update + birthyear + birthmonth +
                       `Mom CMV positive` + `Mom HSV1 positive` + `Mom HSV2 positive` + 
                       `Mom EBV positive` + `Mom Parvo positive` + `Mom HSV6 positive` + 
                       `Mom HSV7 positive` + `Mom RSV positive` + `Mom VZV positive` + `years`
                     , data = curr_model_data, family = "binomial")
  
  curr_summary <- summary(curr_model)
  
  survival_outputs$`Childcare Estimate`[i] <- curr_model$coefficients[2]
  survival_outputs$`Childcare P`[i] <- curr_summary[["coefficients"]][2,4]
  survival_outputs$`Breast Estimate`[i] <- curr_model$coefficients[3]
  survival_outputs$`Breast P`[i] <- curr_summary[["coefficients"]][3,4]
  survival_outputs$`Maternal Age Estimate`[i] <- curr_model$coefficients[4]
  survival_outputs$`Maternal Age P`[i] <- curr_summary[["coefficients"]][4,4]
  survival_outputs$`Preterm Estimate`[i] <- curr_model$coefficients[5]
  survival_outputs$`Preterm P`[i] <- curr_summary[["coefficients"]][5,4]
  survival_outputs$`Parity Estimate`[i] <- curr_model$coefficients[6]
  survival_outputs$`Parity P`[i] <- curr_summary[["coefficients"]][6,4]
  survival_outputs$`Delivery Estimate`[i] <- curr_model$coefficients[7]
  survival_outputs$`Delivery P`[i] <- curr_summary[["coefficients"]][7,4]
  survival_outputs$`Sex Estimate`[i] <- curr_model$coefficients[8]
  survival_outputs$`Sex P`[i] <- curr_summary[["coefficients"]][8,4]
  survival_outputs$`Race 1 Estimate`[i] <- curr_model$coefficients[9]
  survival_outputs$`Race 1 P`[i] <- curr_summary[["coefficients"]][9,4]
  survival_outputs$`Race 3 Estimate`[i] <- curr_model$coefficients[10]
  survival_outputs$`Race 3 P`[i] <- curr_summary[["coefficients"]][10,4]
  survival_outputs$`Race 4 Estimate`[i] <- curr_model$coefficients[11]
  survival_outputs$`Race 4 P`[i] <- curr_summary[["coefficients"]][11,4]
  survival_outputs$`Race 6 Estimate`[i] <- curr_model$coefficients[12]
  survival_outputs$`Race 6 P`[i] <- curr_summary[["coefficients"]][12,4]
  survival_outputs$`Race 7 Estimate`[i] <- curr_model$coefficients[13]
  survival_outputs$`Race 7 P`[i] <- curr_summary[["coefficients"]][13,4]
  survival_outputs$`Race 10 Estimate`[i] <- curr_model$coefficients[14]
  survival_outputs$`Race 10 P`[i] <- curr_summary[["coefficients"]][14,4]
  survival_outputs$`Birth 2001 Estimate`[i] <- curr_model$coefficients[15]
  survival_outputs$`Birth 2001 P`[i] <- curr_summary[["coefficients"]][15,4]
  survival_outputs$`Birth 2002 Estimate`[i] <- curr_model$coefficients[16]
  survival_outputs$`Birth 2002 P`[i] <- curr_summary[["coefficients"]][16,4]
  survival_outputs$`Birth 2003 Estimate`[i] <- curr_model$coefficients[17]
  survival_outputs$`Birth 2003 P`[i] <- curr_summary[["coefficients"]][17,4]
  survival_outputs$`Birth 2004 Estimate`[i] <- curr_model$coefficients[18]
  survival_outputs$`Birth 2004 P`[i] <- curr_summary[["coefficients"]][18,4]
  survival_outputs$`Birth 2005 Estimate`[i] <- curr_model$coefficients[19]
  survival_outputs$`Birth 2005 P`[i] <- curr_summary[["coefficients"]][19,4]
  survival_outputs$`Birth 2006 Estimate`[i] <- curr_model$coefficients[20]
  survival_outputs$`Birth 2006 P`[i] <- curr_summary[["coefficients"]][20,4]
  survival_outputs$`Birth 2008 Estimate`[i] <- curr_model$coefficients[21]
  survival_outputs$`Birth 2008 P`[i] <- curr_summary[["coefficients"]][21,4]
  survival_outputs$`Birth 2009 Estimate`[i] <- curr_model$coefficients[22]
  survival_outputs$`Birth 2009 P`[i] <- curr_summary[["coefficients"]][22,4]
  survival_outputs$`Birth 2010 Estimate`[i] <- curr_model$coefficients[23]
  survival_outputs$`Birth 2010 P`[i] <- curr_summary[["coefficients"]][23,4]
  survival_outputs$`Birth 2011 Estimate`[i] <- curr_model$coefficients[24]
  survival_outputs$`Birth 2011 P`[i] <- curr_summary[["coefficients"]][24,4]
  survival_outputs$`Birth Jan Estimate`[i] <- curr_model$coefficients[25]
  survival_outputs$`Birth Jan P`[i] <- curr_summary[["coefficients"]][25,4]
  survival_outputs$`Birth Feb Estimate`[i] <- curr_model$coefficients[26]
  survival_outputs$`Birth Feb P`[i] <- curr_summary[["coefficients"]][26,4]
  survival_outputs$`Birth Mar Estimate`[i] <- curr_model$coefficients[27]
  survival_outputs$`Birth Mar P`[i] <- curr_summary[["coefficients"]][27,4]
  survival_outputs$`Birth Apr Estimate`[i] <- curr_model$coefficients[28]
  survival_outputs$`Birth Apr P`[i] <- curr_summary[["coefficients"]][28,4]
  survival_outputs$`Birth May Estimate`[i] <- curr_model$coefficients[29]
  survival_outputs$`Birth May P`[i] <- curr_summary[["coefficients"]][29,4]
  survival_outputs$`Birth June Estimate`[i] <- curr_model$coefficients[30]
  survival_outputs$`Birth June P`[i] <- curr_summary[["coefficients"]][30,4]
  survival_outputs$`Birth July Estimate`[i] <- curr_model$coefficients[31]
  survival_outputs$`Birth July P`[i] <- curr_summary[["coefficients"]][31,4]
  survival_outputs$`Birth Aug Estimate`[i] <- curr_model$coefficients[32]
  survival_outputs$`Birth Aug P`[i] <- curr_summary[["coefficients"]][32,4]
  survival_outputs$`Birth Oct Estimate`[i] <- curr_model$coefficients[33]
  survival_outputs$`Birth Oct P`[i] <- curr_summary[["coefficients"]][33,4]
  survival_outputs$`Birth Nov Estimate`[i] <- curr_model$coefficients[34]
  survival_outputs$`Birth Nov P`[i] <- curr_summary[["coefficients"]][34,4]
  survival_outputs$`Birth Dec Estimate`[i] <- curr_model$coefficients[35]
  survival_outputs$`Birth Dec P`[i] <- curr_summary[["coefficients"]][35,4]
  survival_outputs$`Mom CMV Estimate`[i] <- curr_model$coefficients[36]
  survival_outputs$`Mom CMV P`[i] <- curr_summary[["coefficients"]][36,4]
  survival_outputs$`Mom HSV1 Estimate`[i] <- curr_model$coefficients[37]
  survival_outputs$`Mom HSV1 P`[i] <- curr_summary[["coefficients"]][37,4]
  survival_outputs$`Mom HSV2 Estimate`[i] <- curr_model$coefficients[38]
  survival_outputs$`Mom HSV2 P`[i] <- curr_summary[["coefficients"]][38,4]
  survival_outputs$`Mom EBV Estimate`[i] <- curr_model$coefficients[39]
  survival_outputs$`Mom EBV P`[i] <- curr_summary[["coefficients"]][39,4]
  survival_outputs$`Mom Parvo Estimate`[i] <- curr_model$coefficients[40]
  survival_outputs$`Mom Parvo P`[i] <- curr_summary[["coefficients"]][40,4]
  survival_outputs$`Mom HSV6 Estimate`[i] <- curr_model$coefficients[41]
  survival_outputs$`Mom HSV6 P`[i] <- curr_summary[["coefficients"]][41,4]
  survival_outputs$`Mom HSV7 Estimate`[i] <- curr_model$coefficients[42]
  survival_outputs$`Mom HSV7 P`[i] <- curr_summary[["coefficients"]][42,4]
  survival_outputs$`Mom RSV Estimate`[i] <- curr_model$coefficients[43]
  survival_outputs$`Mom RSV P`[i] <- curr_summary[["coefficients"]][43,4]
  survival_outputs$`Mom VZV Estimate`[i] <- curr_model$coefficients[44]
  survival_outputs$`Mom VZV P`[i] <- curr_summary[["coefficients"]][44,4]
  survival_outputs$`Age Estimate`[i] <- curr_model$coefficients[45]
  survival_outputs$`Age P`[i] <- curr_summary[["coefficients"]][45,4]
}

#Calculate BH threshold
BH <- survival_outputs %>% dplyr::select("Virus", ends_with("P")) %>% 
  pivot_longer(cols = contains("P")) %>% mutate(qvalue = p.adjust(value, method = "fdr")) %>% 
  dplyr::select(-value) %>% pivot_wider(names_from = name, values_from = qvalue)

survival_outputs <- survival_outputs %>% dplyr::select(-ends_with("P")) %>%
  left_join(BH)

#volcano plots
childcare <- ggplot() + geom_point(data = survival_outputs, 
                                  aes(x = -(`Childcare Estimate`),
                                      y = -log10(`Childcare P`)),
                                  alpha = .4) +
  theme_light() +
  geom_hline(yintercept = -log10(0.05)) +
  xlab("Coefficient") +
  ylab("-log10(FDR)") +
  labs(title = "Childcare") +
  theme(text = element_text(size = 8))

ggsave("Figures/20250329 - Fig 9 - childcare glm volcano.jpeg",
       childcare,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300)

breast <- ggplot() + geom_point(data = survival_outputs, 
                                  aes(x = (`Breast Estimate`),
                                      y = -log10(`Breast P`)),
                                  alpha = .4) +
  theme_light() +
  geom_hline(yintercept = -log10(0.05)) +
  xlab("Coefficient") +
  ylab("-log10(FDR)") +
  labs(title = "Breastfed") +
  theme(text = element_text(size = 8))

ggsave("Figures/20250329 - Fig 9 - breastfed glm volcano.jpeg",
       breast,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300)

matage <- ggplot() + geom_point(data = survival_outputs, 
                      aes(x = (`Maternal Age Estimate`),
                          y = -log10(`Maternal Age P`)),
                      alpha = .4) +
  theme_light() +
  geom_hline(yintercept = -log10(0.05)) +
  xlab("Coefficient") +
  ylab("-log10(FDR)") +
  labs(title = "Maternal Age") +
  theme(text = element_text(size = 8))

ggsave("Figures/20250329 - Fig 9 - maternal age glm volcano.jpeg",
       matage,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300)

preterm <- ggplot() + geom_point(data = survival_outputs, 
                      aes(x = (`Preterm Estimate`),
                          y = -log10(`Preterm P`)),
                      alpha = .4) +
  theme_light() +
  geom_hline(yintercept = -log10(0.05)) +
  xlab("Coefficient") +
  ylab("-log10(FDR)") +
  labs(title = "Preterm") +
  theme(text = element_text(size = 8))

ggsave("Figures/20250329 - Fig 9 - preterm glm volcano.jpeg",
       preterm,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300)

Parity <- ggplot() + geom_point(data = survival_outputs, 
                      aes(x = (`Parity Estimate`),
                          y = -log10(`Parity P`)),
                      alpha = .4) +
  theme_light() +
  geom_hline(yintercept = -log10(0.05)) +
  xlab("Coefficient") +
  ylab("-log10(FDR)") +
  labs(title = "Maternal Parity > 0") +
  theme(text = element_text(size = 8))

ggsave("Figures/20250329 - Fig 9 - parity glm volcano.jpeg",
       Parity,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300)

Delivery <- ggplot() + geom_point(data = survival_outputs, 
                      aes(x = (`Delivery Estimate`),
                          y = -log10(`Delivery P`)),
                      alpha = .4) +
  theme_light() +
  geom_hline(yintercept = -log10(0.05)) +
  xlab("Coefficient") +
  ylab("-log10(FDR)") +
  labs(title = "Delivery") +
  theme(text = element_text(size = 8))

ggsave("Figures/20250329 - Fig 9 - delivery glm volcano.jpeg",
       Delivery,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300)

Sex <- ggplot() + geom_point(data = survival_outputs, 
                      aes(x = (`Sex Estimate`),
                          y = -log10(`Sex P`)),
                      alpha = .4) +
  theme_light() +
  geom_hline(yintercept = -log10(0.05)) +
  xlab("Coefficient") +
  ylab("-log10(FDR)") +
  labs(title = "Sex") +
  theme(text = element_text(size = 8))

ggsave("Figures/20250329 - Fig 9 - sex glm volcano.jpeg",
       Sex,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300)

Race <- ggplot() + geom_point(data = survival_outputs, 
                      aes(x = (`Race 1 Estimate`),
                          y = -log10(`Race 1 P`)),
                      alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Race 3 Estimate`),
                 y = -log10(`Race 3 P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Race 4 Estimate`),
                 y = -log10(`Race 4 P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Race 6 Estimate`),
                 y = -log10(`Race 6 P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Race 7 Estimate`),
                 y = -log10(`Race 7 P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Race 10 Estimate`),
                 y = -log10(`Race 10 P`)),
             alpha = .2) +
  theme_light() +
  geom_hline(yintercept = -log10(0.05)) +
  xlab("Coefficient") +
  ylab("-log10(FDR)") +
  labs(title = "Race") +
  theme(text = element_text(size = 8))

ggsave("Figures/20250329 - Fig 9 - race glm volcano.jpeg",
       Race,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300)

Birth_Year <- ggplot() + geom_point(data = survival_outputs, 
                              aes(x = (`Birth 2001 Estimate`),
                                  y = -log10(`Birth 2001 P`)),
                              alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth 2002 Estimate`),
                 y = -log10(`Birth 2002 P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth 2003 Estimate`),
                 y = -log10(`Birth 2003 P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth 2004 Estimate`),
                 y = -log10(`Birth 2004 P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth 2005 Estimate`),
                 y = -log10(`Birth 2005 P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth 2006 Estimate`),
                 y = -log10(`Birth 2006 P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth 2008 Estimate`),
                 y = -log10(`Birth 2008 P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth 2009 Estimate`),
                 y = -log10(`Birth 2009 P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth 2010 Estimate`),
                 y = -log10(`Birth 2010 P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth 2011 Estimate`),
                 y = -log10(`Birth 2011 P`)),
             alpha = .2) +
  theme_light() +
  geom_hline(yintercept = -log10(0.05)) +
  xlab("Coefficient") +
  ylab("-log10(FDR)") +
  labs(title = "Birth Year") +
  theme(text = element_text(size = 8))

ggsave("Figures/20250329 - Fig 9 - year glm volcano.jpeg",
       Birth_Year,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300)


Birth_Month <- ggplot() + geom_point(data = survival_outputs, 
                                    aes(x = (`Birth Jan Estimate`),
                                        y = -log10(`Birth Jan P`)),
                                    alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth Feb Estimate`),
                 y = -log10(`Birth Feb P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth Mar Estimate`),
                 y = -log10(`Birth Mar P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth Apr Estimate`),
                 y = -log10(`Birth Apr P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth May Estimate`),
                 y = -log10(`Birth May P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth June Estimate`),
                 y = -log10(`Birth June P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth July Estimate`),
                 y = -log10(`Birth July P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth Aug Estimate`),
                 y = -log10(`Birth Aug P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth Oct Estimate`),
                 y = -log10(`Birth Oct P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth Nov Estimate`),
                 y = -log10(`Birth Nov P`)),
             alpha = .2) +
  geom_point(data = survival_outputs, 
             aes(x = (`Birth Dec Estimate`),
                 y = -log10(`Birth Dec P`)),
             alpha = .2) +
  theme_light() +
  geom_hline(yintercept = -log10(0.05)) +
  xlab("Coefficient") +
  ylab("-log10(FDR)") +
  labs(title = "Birth Month") +
  theme(text = element_text(size = 8))

ggsave("Figures/20250329 - Fig 9 - month glm volcano.jpeg",
       Birth_Month,
       scale = 1,
       height = 2,
       width = 2,
       units = "in",
       dpi = 300)

ggplot() + geom_point(data = survival_outputs, 
                      aes(x = (`Age Estimate`),
                          y = -log10(`Age P`)),
                      alpha = .4) +
  theme_light() +
  geom_hline(yintercept = -log10(0.05)) +
  xlab("Coefficient") +
  ylab("-log10(FDR)") +
  labs(title = "Age") +
  theme(text = element_text(size = 8))
