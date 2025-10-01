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

# a single serotype of each species (Best represented)
hsv_VAR <- BBC_VARscores %>%
  filter(taxon_species %in% c("Human adenovirus A serotype 12", 
                             "Human adenovirus B serotype 16",
                              "Human adenovirus C serotype 2",
                              "Human adenovirus D37",
                              "Human adenovirus E serotype 4", 
                              "Human adenovirus F serotype 40")) %>%
 mutate(taxon_species = as.factor(taxon_species)) %>%
  arrange(taxon_species)

dotplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>% group_by(Plasma_type),
             aes(x = interaction(Plasma_type, taxon_species), y = vir_score,
                 color = Plasma_type),
             position = position_quasirandom(),
             size = .1,
             alpha = .3) +
  ylab("VARscore") +
  xlab("Virus") +
  geom_hline(yintercept = 1.936591) +
  scale_x_discrete(labels=c("cord.Human adenovirus A serotype 12" = "Adeno A",
                            "child.Human adenovirus A serotype 12" = "Adeno A",
                            "cord.Human adenovirus B serotype 16" = "Adeno B",
                            "child.Human adenovirus B serotype 16" = "Adeno B",
                            "cord.Human adenovirus C serotype 2" = "Adeno C",
                            "child.Human adenovirus C serotype 2" = "Adeno C",
                            "cord.Human adenovirus D37" = "Adeno D",
                            "child.Human adenovirus D37" = "Adeno D",
                            "cord.Human adenovirus E serotype 4" = "Adeno E",
                            "child.Human adenovirus E serotype 4" = "Adeno E",
                            "cord.Human adenovirus F serotype 40" = "Adeno F",
                            "child.Human adenovirus F serotype 40" = "Adeno F"
                            )) +
  theme_light() +
  scale_color_manual(values = c("black", "red")) +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - adenovirus VARscore dotplot.jpeg",
       dotplot,
       scale = 1,
       height = 2,
       width = 7.5,
       units = "in",
       dpi = 300)

timeplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% filter(Plasma_type == "child") %>%
               mutate(taxon_species = case_when(
                 taxon_species == "Human adenovirus A serotype 12" ~ "Adeno A",
                 taxon_species == "Human adenovirus B serotype 16" ~ "Adeno B",
                 taxon_species == "Human adenovirus C serotype 2" ~ "Adeno C",
                 taxon_species == "Human adenovirus D37" ~ "Adeno D",
                 taxon_species == "Human adenovirus E serotype 4" ~ "Adeno E",
                 taxon_species == "Human adenovirus F serotype 40" ~ "Adeno F",
               )),
             aes(x = years, y = vir_score),
             size = .3,
             alpha = .3) +
  ylab("VARscore") +
  xlab("Age (years)") +
  theme_light() + 
  geom_hline(yintercept = 1.936591) +
  facet_wrap(vars(taxon_species)) + 
  theme(text = element_text(size = 8))

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - child adenovirus VARscore vs time dotplot.jpeg",
       timeplot,
       scale = 1,
       height = 3,
       width = 4,
       units = "in",
       dpi = 300)

ggplot() + 
  geom_point(data = hsv_VAR %>% filter(Plasma_type == "child"),
             aes(x = `Plasma_Date`, y = vir_score),
             size = 1,
             alpha = .3) +
  ylab("VARscore") +
  xlab("Date") +
  theme_light() + 
  geom_smooth(data = hsv_VAR %>% filter(Plasma_type == "child"),
              aes(x = `Plasma_Date`, y = vir_score),
              formula = y ~ x,
              span = 0.25) +
  facet_wrap(vars(taxon_species))

hsv_serobyage = hsv_VAR %>% filter(Plasma_type == "child") %>%
  mutate(age_group = case_when(
    years < 1 ~ "0 - 1",
    years < 2 ~ "1 - 2",
    years < 3 ~ "2 - 3",
    years < 4 ~ "3 - 4",
    years >= 4 ~ "4+"
  )) %>%
  mutate(seropositive = (vir_score > 1.936591)) %>%
  ungroup() %>% group_by(age_group, taxon_species) %>%
  summarise(proportion = sum(seropositive)/n(), number = n())
sero_byage <- ggplot() + 
  geom_col(data = hsv_serobyage,
           aes(x = taxon_species, y = proportion,
               fill = age_group),
           position = position_dodge(),
           color = "black"
  ) +
  ylab("Proportion positive") +
  scale_x_discrete(labels=c("Human adenovirus A serotype 12" = "Adeno A",
                            "Human adenovirus B serotype 16" = "Adeno B",
                            "Human adenovirus C serotype 2" = "Adeno C",
                            "Human adenovirus D37" = "Adeno D",
                            "Human adenovirus E serotype 4" = "Adeno E",
                            "Human adenovirus F serotype 40" = "Adeno F"
  )) +
  theme_light() + 
  scale_fill_viridis_d() +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - child adenosvirus seropositivity by age.jpeg",
       sero_byage,
       scale = 1,
       height = 1.5,
       width = 3.5,
       units = "in",
       dpi = 300)

##############
# cross reactivity in kids
##############

HSV1_HSV2 <- hsv_VAR %>% ungroup() %>%
  pivot_wider(id_cols = c(ID, Plasma_type), names_from = taxon_species,
              values_from = vir_score, values_fn = max)

HSV_xreact <- ggplot() + 
  geom_point(data = HSV1_HSV2,
             aes(x = `Human adenovirus B serotype 16`, 
                 y = `Human adenovirus C serotype 2`,
                 color = Plasma_type),
             size = .3,
             alpha = .3) +
  xlab("Adeno B") +
  ylab("Adeno C") +
  theme_light() +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(vars(Plasma_type)) +
  theme(text = element_text(size = 8),
        legend.position = "none")

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - adeno B C cross reactivity.jpeg",
       HSV_xreact,
       scale = 1,
       height = 1.5,
       width = 3,
       units = "in",
       dpi = 300)

HSV1_HSV2 <- HSV1_HSV2 %>% ungroup() %>%
  filter(Plasma_type == "child")

library(GGally)

pairs_dots <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
    geom_point(..., alpha = 0.2) 
}

ggpairs(data = HSV1_HSV2,
        columns = c(3:8),
        #mapping = aes(alpha = .1)
        lower = list(continuous = pairs_dots)) +
  theme_light()
