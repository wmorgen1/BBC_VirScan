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
########
########
# Data Prep
#######
library(readxl)
sample_meta <- read_excel("BBC/sample_meta.xls") %>% filter(ID %nin% c(4452, 8754))

viruses_totest <- read.csv("BBC/Data/20240424 - BBC cord child volcano.csv")

sample_meta <- sample_meta %>% mutate(age = difftime(time1 = Plasma_Date, time2 = DOD, units = "days")) %>%
  mutate(years = (age %>% as.numeric())/365.25) %>% mutate(Plasma_type = ifelse(Plasma_type=="Postnatal","child",Plasma_type))

breast_fed_meta <- read_excel("BBC/breastfed_update ID with codebook.xlsx")
parity_meta <- read_excel("BBC/parity.xlsx")
care_meta <- read_csv("BBC/child care.csv") %>% filter((ID != 1574 | I_childcare == 1)) %>%
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
additional_mom_meta <- read_excel("BBC/pheno_update2.xlsx")

sample_meta <- sample_meta %>% left_join(breast_fed_meta) %>%
  left_join(parity_meta) %>% left_join(care_meta) %>%
  left_join(additional_mom_meta)

BBC_VARscores <- read_rds("BBC/Data/20240424 - BBC VARscores - maternal abs removed.rds")

BBC_VARscores <- BBC_VARscores %>% dplyr::rename(vir_score = value) %>%
  left_join(sample_meta, multiple = "all")

hsv_VAR <- BBC_VARscores %>%
  filter(taxon_species == "Epstein-Barr virus") %>%
  mutate(race_update = as.factor(race_update)) %>%
  #mutate(race_update = fct_relevel(race_update, c("2", "1", "3", "4", "5", "6", "7", "10")))%>%
  mutate(birthyear = as.factor(year(DOD))) %>%
  mutate(birthyear = fct_relevel(birthyear, c("2007", "2000", "2001", "2002",
                                              "2003", "2004", "2005", "2006",
                                              "2008", "2009", "2010", "2011")))%>%
  mutate(birthmonth = as.factor(month(DOD))) %>%
  mutate(birthmonth = fct_relevel(birthmonth, c("9", "1", "2", "3",
                                                "4", "5", "6", "7",
                                                "8", "10", "11", "12"))) %>%
  mutate(race_update = ifelse(race_update == "5", "1",race_update))

child_meta <- hsv_VAR %>% filter(Plasma_type == "child") 

library(table1)

table1(~ years + mage + factor(birthmonth) + factor(birthyear) + factor(sex) + factor(race_update) + factor(parity2) + 
         factor(DELTYPE) + factor(preterm1) + factor(breast_fed1) + geaa_o, data = child_meta)

## Codebook
#years - age of child at time of plasma collection in years
#mage - maternal age at childbirth
#sex - 0 = F, 1 = M
#Race - 1 = Black/African American, 2 = White, 3 = Hispanic, 4 = Asian, 6 = Cape Verdian, 7 = Pacific Islander, 8 = Mixed Race
#Parity2 - 0 = nulliparity, 1 = parity > 0
#DELTYPE - 1 = vaginal, 1 = C-section
#preterm1 - 0 = term, 1 = preterm
#breast_fed1 - 0 = bottle-fed, 1 = both, 2 = exclusively breastfed
#birth month - 1 = Jan, 2 = Feb ... etc
#birth year
