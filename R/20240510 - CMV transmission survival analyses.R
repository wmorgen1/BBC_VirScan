setwd("C:/Users/wmorgen1/Documents/Larman_Lab/")
library(tidyverse)
library(RColorBrewer)
library(gtable)
library(grid)
library(ggplot2)
library(ggbeeswarm)
library(data.table)
#library(ARscore)

########
# Data Prep
#######
library(readxl)
sample_meta <- read_excel("BBC/sample_meta.xls")

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

#####
# Regression and Survival Analyses
#####

CMV_model_data <- hsv_VAR %>% ungroup() %>% filter(Plasma_type == "child") %>%
  left_join(Maternal_pos) %>% filter(taxon_species == "Human cytomegalovirus") %>%
  mutate(`Mom CMV positive` = as.factor(`Mom CMV positive`)) %>%
  #filter(`Mom CMV positive` == TRUE) %>% filter(breast_fed1 != "NA") %>%
  mutate(`kid pos` = as.factor(vir_score > 1.936591)) %>%
  mutate(cc = case_when(
    (years < begin_of_cc) | (I_childcare == 0) ~ "No Childcare",
    years > begin_of_cc ~ "Childcare",
    T ~ NA
  )) %>%
  mutate(breast_fed = as.factor(breast_fed1 %in% c(1, 2))) 
  ## investigate factors

object <- glm(`kid pos` ~ `Mom CMV positive` + `cc` + `breast_fed` + mage + preterm1 + 
                parity2 + DELTYPE + sex + race_update, data = CMV_model_data, family = "binomial")
summary(object) # results in logit coefficients
exp(cbind(OR = coef(object), confint(object))) # results in ORs and their CIs

  ## only significant factors
CMV_model_data <- CMV_model_data %>% mutate(`Mom CMV positive` = as.factor(`Mom CMV positive`)) %>%
  mutate(`Mom CMV positive` = fct_relevel(`Mom CMV positive`, c('FALSE', 'TRUE'))) %>% 
  mutate(`Breast fed` = as.factor(`breast_fed`)) %>%
  mutate(`Breast fed` = fct_relevel(`Breast fed`, c('FALSE', 'TRUE')))

  
object <- glm(`kid pos` ~ `cc` + `breast_fed`:`Mom CMV positive`, 
              data = CMV_model_data, family = "binomial")
summary(object) # results in logit coefficients
exp(cbind(OR = coef(object), confint(object))) # results in ORs and their CIs


ORs <- sjPlot::plot_model(object,
                   axis.lim = c(0.05, 2)) +
  theme_light() +
  scale_y_log10(limits = c(0.025, 2),
                     breaks = c(0.033, 0.1, 0.33, 1),
                minor_breaks = NULL) +
  theme(plot.title = element_blank(),
        text = element_text(size = 8)) +
  geom_hline(yintercept = 1)

#ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240515 - Logistic Regression ORs.jpeg",
#       ORs,
#       scale = 1,
#       height = 1.5,
#       width = 4,
#       units = "in",
#       dpi = 300)

library(survival)
library(survminer)
library(icenReg)
library(ggplot2)

CMV_model_data <- CMV_model_data %>% mutate(years_0 = ifelse(`kid pos` == TRUE, 0, years)) %>%
  mutate(years_1 = ifelse(`kid pos` == TRUE, years, Inf))

CMV_model_data_2 <- CMV_model_data %>% mutate(years_0 = ifelse(`kid pos` == TRUE, 0, years)) %>%
  mutate(years_1 = ifelse(`kid pos` == TRUE, years, Inf)) %>% mutate(feed_and_mom = interaction(`Mom CMV positive`, `breast_fed`))

#res.cox <- ic_sp(Surv(time = years_0, time2 = years_1, type = 'interval2') ~ `feed_and_mom`,
#                 model = "ph",
#                 data = CMV_model_data_2 %>% filter(!is.na(`Mom CMV positive`)) %>%
#                   #filter(!is.na(`cc`)) %>%
#                   filter(!is.na(`breast_fed`)),
#                 bs_samples = 500) 

#summary(res.cox)
#plot(res.cox)
#newdata <- data.frame('feed_and_mom' = c('TRUE.TRUE', 'TRUE.FALSE', 'FALSE.TRUE', 'FALSE.FALSE'),
#                      'cc' = c("No Childcare", "No Childcare", "No Childcare", "No Childcare"))

#names(newdata) = c('feed_and_mom', 'cc')
#newdata <- newdata %>% mutate(feed_and_mom = as.factor(feed_and_mom))
#plot(res.cox, newdata)

library(interval)

ictest(Surv(time = years_0, time2 = years_1, type = 'interval2') ~ `feed_and_mom`,
       data = CMV_model_data_2)

ictest(Surv(time = years_0, time2 = years_1, type = 'interval2') ~ `cc`,
       data = CMV_model_data_2)


ggsurvplot(survfit(Surv(time = years_0, time2 = years_1, type = 'interval2') ~ `feed_and_mom`,
                   data = CMV_model_data_2), data = CMV_model_data_2,
           size = 1,                 # change line size
           linetype = c(2, 2, 2, 1),
           palette =
             c("black", "#E7B800", "#2E9FDF", "red"),# custom color palettes
          legend.labs =
             c("Mom-, BF-", "Mom+, BF-", "Mom-, BF+", "Mom+, BF+"),    # Change legend labels
          risk.table = F,        # Add risk table
          risk.table.col = "strata",# Risk table color by groups
           risk.table.height = 0.25, # Useful to change when you have multiple groups
           ggtheme = theme_bw(),
          font.main = c(10, "plain", "black"),
          font.x = c(8, "plain", "black"),
          font.y = c(8, "plain", "black"),
          font.tickslab = c(8, "plain", "black"),
          risk.table.fontsize = 3)


ggsave("BBC/Figures/20241110 - CMV Survival Curve - big.jpeg",
       scale = 1,
       height = 6,
       width = 6,
       units = "in",
       dpi = 300)
  

ggsurvplot(survfit(Surv(time = years_0, time2 = years_1, type = 'interval2') ~ `cc`,
                   data = CMV_model_data_2 %>% filter(feed_and_mom != "TRUE.TRUE")), 
           data = CMV_model_data_2 %>% filter(feed_and_mom != "TRUE.TRUE"),
           size = 1,                 # change line size
           #palette =
           #   c("#E7B800", "#2E9FDF", "red", "black"),# custom color palettes
           #legend.labs =
          #    c("1", "2", "3", "4"),    # Change legend labels
           risk.table.height = 0.25, # Useful to change when you have multiple groups
           ggtheme = theme_bw())