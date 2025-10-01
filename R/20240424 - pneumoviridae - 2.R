setwd("C:/Users/wmorgen1/Documents/Larman_Lab/BBC/")
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
sample_meta <- read_excel("sample_meta.xls")

sample_meta <- sample_meta %>% mutate(age = difftime(time1 = Plasma_Date, time2 = DOD, units = "days")) %>%
  mutate(years = (age %>% as.numeric())/365.25) %>% mutate(Plasma_type = ifelse(Plasma_type=="Postnatal","child",Plasma_type))

BBC_VARscores <- read_rds("Data/20240424 - BBC VARscores - maternal abs removed.rds")

BBC_VARscores <- BBC_VARscores %>% dplyr::rename(vir_score = value) %>%
  left_join(sample_meta, multiple = "all")

# excluding unspecified
hsv_VAR <- BBC_VARscores %>%
  filter(taxon_genus %in% c("Orthopneumovirus", "Metapneumovirus")) %>%
  filter(taxon_species != "Human respiratory syncytial virus") %>%
  mutate(taxon_species = as.factor(taxon_species)) %>%
  arrange(taxon_species)

hsv_VAR <- hsv_VAR %>% mutate(time_of_year = case_when(
  month(Plasma_Date) == 1 ~ 0 + day(Plasma_Date)/31,
  month(Plasma_Date) == 2 ~ 1 + day(Plasma_Date)/28,
  month(Plasma_Date) == 3 ~ 2 + day(Plasma_Date)/31,
  month(Plasma_Date) == 4 ~ 3 + day(Plasma_Date)/30,
  month(Plasma_Date) == 5 ~ 4 + day(Plasma_Date)/31,
  month(Plasma_Date) == 6 ~ 5 + day(Plasma_Date)/30,
  month(Plasma_Date) == 7 ~ 6 + day(Plasma_Date)/31,
  month(Plasma_Date) == 8 ~ 7 + day(Plasma_Date)/31,
  month(Plasma_Date) == 9 ~ 8 + day(Plasma_Date)/30,
  month(Plasma_Date) == 10 ~ 9 + day(Plasma_Date)/31,
  month(Plasma_Date) == 11 ~ 10 + day(Plasma_Date)/30,
  month(Plasma_Date) == 12 ~ 11 + day(Plasma_Date)/31
))

dateplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% filter(Plasma_type == "child")%>%
               mutate(taxon_species = case_when(
                 taxon_species == "Human metapneumovirus" ~ "MPnV",
                 taxon_species == "Human respiratory syncytial virus A" ~ "RSV A",
                 taxon_species == "Human respiratory syncytial virus B" ~ "RSV B"
               )),
             aes(x = `time_of_year`, y = vir_score),
             size = 1,
             alpha = .3) +
  ylab("VARscore") +
  xlab("Month") +
  theme_light() + 
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) +
  geom_smooth(data = hsv_VAR %>% filter(Plasma_type == "child")%>%
                mutate(taxon_species = case_when(
                  taxon_species == "Human metapneumovirus" ~ "MPnV",
                  taxon_species == "Human respiratory syncytial virus A" ~ "RSV A",
                  taxon_species == "Human respiratory syncytial virus B" ~ "RSV B"
                )),
              aes(x = `time_of_year`, y = vir_score),
              formula = y ~ x,
              span = 0.5) +
  facet_wrap(vars(taxon_species))

ggsave("Figures/20241110 - child pneumo VARscore vs month dotplot - seasonality and outbreaks.jpeg",
       dateplot,
       scale = 1,
       height = 2,
       width = 4,
       units = "in",
       dpi = 300)


dateplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% filter(Plasma_type == "cord")%>%
               mutate(taxon_species = case_when(
                 taxon_species == "Human metapneumovirus" ~ "MPnV",
                 taxon_species == "Human respiratory syncytial virus A" ~ "RSV A",
                 taxon_species == "Human respiratory syncytial virus B" ~ "RSV B"
               )),
             aes(x = `time_of_year`, y = vir_score),
             size = 1,
             alpha = .3) +
  ylab("VARscore") +
  xlab("Month") +
  theme_light() + 
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) +
  geom_smooth(data = hsv_VAR %>% filter(Plasma_type == "cord")%>%
                mutate(taxon_species = case_when(
                  taxon_species == "Human metapneumovirus" ~ "MPnV",
                  taxon_species == "Human respiratory syncytial virus A" ~ "RSV A",
                  taxon_species == "Human respiratory syncytial virus B" ~ "RSV B"
                )),
              aes(x = `time_of_year`, y = vir_score),
              formula = y ~ x,
              span = 0.5) +
  facet_wrap(vars(taxon_species))

ggsave("Figures/20241110 - cord pneumo VARscore vs month dotplot - seasonality and outbreaks.jpeg",
       dateplot,
       scale = 1,
       height = 2,
       width = 4,
       units = "in",
       dpi = 300)

