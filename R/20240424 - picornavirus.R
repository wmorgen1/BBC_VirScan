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

BBC_VARscores <- read_rds("BBC/Data/20240424 - BBC VARscores - maternal abs removed.rds")

BBC_VARscores <- BBC_VARscores %>% dplyr::rename(vir_score = value) %>%
  left_join(sample_meta, multiple = "all")

hsv_VAR <- BBC_VARscores %>%
  filter(taxon_genus %in% c("Enterovirus", "Cosavirus", "Kobuvirus", "Parechovirus", "Salivivirus")) %>%
  filter(taxon_species %nin% c("Canine kobuvirus US-PC0082")) %>% # not canine viruses
  mutate(taxon_species = as.factor(taxon_species)) %>%
  mutate(taxon_species = fct_relevel(taxon_species, "Enterovirus A", "Enterovirus B",
                                     "Enterovirus C", "Enterovirus D", 
                                     "Rhinovirus A", "Rhinovirus B", "Human parechovirus 1",
                                     "Human parechovirus 2", "Cosavirus A", "Aichi virus")) %>%
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
  scale_x_discrete(labels=c("cord.Enterovirus A" = "EVA",
                            "child.Enterovirus A" = "EVA",
                            "cord.Enterovirus B" = "EVB",
                            "child.Enterovirus B" = "EVB",
                            "cord.Enterovirus C" = "EVC",
                            "child.Enterovirus C" = "EVC",
                            "cord.Enterovirus D" = "EVD",
                            "child.Enterovirus D" = "EVD",
                            "cord.Rhinovirus A" = "RVA",
                            "child.Rhinovirus A" = "RVA",
                            "cord.Rhinovirus B" = "RVB",
                            "child.Rhinovirus B" = "RVB",
                            "cord.Human parechovirus 1" = "pEV1",
                            "child.Human parechovirus 1" = "pEV1",
                            "cord.Human parechovirus 2" = "pEV2",
                            "child.Human parechovirus 2" = "pEV2",
                            "cord.Cosavirus A" = "CVA",
                            "child.Cosavirus A" = "CVA",
                            "cord.Aichi virus" = "Aichi",
                            "child.Aichi virus" = "Aichi"
                            )) +
  theme_light() +
scale_color_manual(values = c("black", "red")) +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - picornavirus VARscore dotplot.jpeg",
       dotplot,
       scale = 1,
       height = 2,
       width = 7.5,
       units = "in",
       dpi = 300)

timeplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% ungroup() %>%
               filter(Plasma_type == "child") %>%
               filter(taxon_species %nin% c("Cosavirus A", "Aichi virus")) %>%
               mutate(taxon_species = case_when(
                 taxon_species == "Enterovirus A" ~ "EVA",
                 taxon_species == "Enterovirus B" ~ "EVB",
                 taxon_species == "Enterovirus C" ~ "EVC",
                 taxon_species == "Enterovirus D" ~ "EVD",
                 taxon_species == "Rhinovirus A" ~ "RVA",
                 taxon_species == "Rhinovirus B" ~ "RVB",
                 taxon_species == "Human parechovirus 1" ~ "pEV1",
                 taxon_species == "Human parechovirus 2" ~ "pEV2"
               )),
             aes(x = years, y = vir_score),
             size = .3,
             alpha = .3) +
  geom_hline(yintercept = 1.936591) +
  ylab("VARscore") +
  xlab("Age (years)") +
  theme_light() + 
  facet_wrap(vars(taxon_species), nrow = 2) + 
  theme(text = element_text(size = 8))

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - child picornavirus VARscore vs time dotplot.jpeg",
       timeplot,
       scale = 1,
       height = 3,
       width = 5,
       units = "in",
       dpi = 300)

dateplot <- ggplot() + 
  geom_point(data = hsv_VAR %>% filter(Plasma_type == "child")  %>%
               filter(taxon_species %nin% c("Cosavirus A", "Aichi virus")) %>%
               mutate(taxon_species = case_when(
                 taxon_species == "Enterovirus A" ~ "EVA",
                 taxon_species == "Enterovirus B" ~ "EVB",
                 taxon_species == "Enterovirus C" ~ "EVC",
                 taxon_species == "Enterovirus D" ~ "EVD",
                 taxon_species == "Rhinovirus A" ~ "RVA",
                 taxon_species == "Rhinovirus B" ~ "RVB",
                 taxon_species == "Human parechovirus 1" ~ "pEV1",
                 taxon_species == "Human parechovirus 2" ~ "pEV2"
               )),
             aes(x = `Plasma_Date`, y = vir_score),
             size = 1,
             alpha = .3) +
  ylab("VARscore") +
  xlab("Date") +
  theme_light() + 
  geom_smooth(data = hsv_VAR %>% filter(Plasma_type == "child")  %>%
                filter(taxon_species %nin% c("Cosavirus A", "Aichi virus")) %>%
                mutate(taxon_species = case_when(
                  taxon_species == "Enterovirus A" ~ "EVA",
                  taxon_species == "Enterovirus B" ~ "EVB",
                  taxon_species == "Enterovirus C" ~ "EVC",
                  taxon_species == "Enterovirus D" ~ "EVD",
                  taxon_species == "Rhinovirus A" ~ "RVA",
                  taxon_species == "Rhinovirus B" ~ "RVB",
                  taxon_species == "Human parechovirus 1" ~ "pEV1",
                  taxon_species == "Human parechovirus 2" ~ "pEV2"
                )),
              aes(x = `Plasma_Date`, y = vir_score),
              formula = y ~ x,
              method = "loess",
              span = 0.25) +
  facet_wrap(vars(taxon_species), nrow = 2)

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - child picorna VARscore vs date dotplot - seasonality and outbreaks.jpeg",
       dateplot,
       scale = 1,
       height = 3,
       width = 7,
       units = "in",
       dpi = 300)

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
  geom_col(data = hsv_serobyage %>% 
               filter(taxon_species %in% c("Enterovirus A", "Enterovirus B",
                                           "Enterovirus C", "Enterovirus D", 
                                           "Rhinovirus A", "Rhinovirus B", "Human parechovirus 1",
                                           "Human parechovirus 2")),
           aes(x = taxon_species, y = proportion,
               fill = age_group),
           position = position_dodge(),
             color = "black"
             ) +
  scale_x_discrete(labels=c("Enterovirus A" = "EVA",
                            "Enterovirus B" = "EVB",
                            "Enterovirus C" = "EVC",
                            "Enterovirus D" = "EVD",
                            "Rhinovirus A" = "RVA",
                            "Rhinovirus B" = "RVB",
                            "Human parechovirus 1" = "pEV1",
                            "Human parechovirus 2" = "pEV2"
  )) +
  ylab("Proportion positive") +
  theme_light() + 
  scale_fill_viridis_d() +
  theme(text = element_text(size = 8),
        axis.title.x = element_blank())

ggsave("C:/Users/wmorgen1/Documents/R/Larman_Lab/BBC/Figures/20240424 - child picornavirus seropositivity by age.jpeg",
       sero_byage,
       scale = 1,
       height = 1.5,
       width = 5,
       units = "in",
       dpi = 300)

##############
# cross reactivity in kids
##############

HSV1_HSV2 <- hsv_VAR %>% ungroup() %>%
  filter(taxon_genus == "Enterovirus") %>%
  pivot_wider(id_cols = c(ID, Plasma_type), names_from = taxon_species,
              values_from = vir_score, values_fn = max)

xreact <- ggplot() + 
  geom_point(data = HSV1_HSV2,
             aes(x = `Enterovirus A`, 
                 y = `Enterovirus B`,
                 color = Plasma_type),
             size = .3,
             alpha = .3) +
  xlab("EVA") +
  ylab("EVB") +
  theme_light() +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(vars(Plasma_type), ncol = 2) +
  theme(text = element_text(size = 8),
        legend.position = "none")

ggsave("BBC/Figures/20240424 - EVA EVB cross reactivity.jpeg",
       xreact,
       scale = 1,
       height = 1.5,
       width = 2.5,
       units = "in",
       dpi = 300)

xreact_2 <- ggplot() + 
  geom_point(data = HSV1_HSV2,
             aes(x = `Enterovirus A`, 
                 y = `Enterovirus B`,
                 color = Plasma_type),
             size = .3,
             alpha = .3) +
  xlab("EVA") +
  ylab("EVB") +
  theme_light() +
  scale_color_manual(values = c("red", "black")) +
  facet_wrap(vars(Plasma_type), ncol = 1) +
  theme(text = element_text(size = 8),
        legend.position = "none")

library(GGally)

HSV1_HSV2 <- HSV1_HSV2 %>% ungroup() %>%
  filter(Plasma_type == "child")

#names(HSV1_HSV2) <- c("ID", "Plasma_type", 
#                      "HSV1", "HSV2",
#                      "VZV", "EBV",
#                      "CMV", "HSV6A", "HSV6B",
#                      "HSV7", "HSV8")

pairs_dots <- function(data, mapping, ...) {
  ggplot(data = data, mapping=mapping) +
    geom_point(..., size = .1, alpha = 0.2) 
}

ggpairs(data = HSV1_HSV2,
        columns = c(3:8),
        #mapping = aes(alpha = .1)
        lower = list(continuous = pairs_dots)) +
  theme_light() + 
  theme(text = element_text(size = 8))

EVA <- lm(`Enterovirus A` ~ `Enterovirus B` + `Enterovirus C` + `Enterovirus D`  + `Rhinovirus A` + `Rhinovirus B` - 1, data = HSV1_HSV2)
EVB <- lm(`Enterovirus B` ~ `Enterovirus A` + `Enterovirus C` + `Enterovirus D`  + `Rhinovirus A` + `Rhinovirus B` - 1, data = HSV1_HSV2)
EVC <- lm(`Enterovirus C` ~ `Enterovirus B` + `Enterovirus A` + `Enterovirus D`  + `Rhinovirus A` + `Rhinovirus B` - 1, data = HSV1_HSV2)
EVD <- lm(`Enterovirus D` ~ `Enterovirus B` + `Enterovirus C` + `Enterovirus A`  + `Rhinovirus A` + `Rhinovirus B` - 1, data = HSV1_HSV2)
RVA <- lm(`Rhinovirus A` ~ `Enterovirus B` + `Enterovirus C` + `Enterovirus D`  + `Enterovirus A` + `Rhinovirus B` - 1, data = HSV1_HSV2)
RVB <- lm(`Rhinovirus B` ~ `Enterovirus B` + `Enterovirus C` + `Enterovirus D`  + `Rhinovirus A` + `Enterovirus A` - 1, data = HSV1_HSV2)

HSV1_HSV2_deconv <- HSV1_HSV2 %>% 
  mutate(`EVA adjusted` = EVA$residuals + `Enterovirus A`) %>% 
  mutate(`EVB adjusted` = EVB$residuals + `Enterovirus B`) %>% 
  mutate(`EVC adjusted` = EVC$residuals + `Enterovirus C`) %>% 
  mutate(`EVD adjusted` = EVD$residuals + `Enterovirus D`) %>% 
  mutate(`RVA adjusted` = RVA$residuals + `Rhinovirus A`) %>% 
  mutate(`RVB adjusted` = RVB$residuals + `Rhinovirus B`)

ggpairs(data = HSV1_HSV2_deconv,
        columns = c(9:14),
        #mapping = aes(alpha = .1)
        lower = list(continuous = pairs_dots)) +
  theme_light() + 
  theme(text = element_text(size = 8))

HSV1_HSV2_deconv <- HSV1_HSV2 %>% 
  mutate(`EVA adjusted` = `Enterovirus A` - (0.789 * `Enterovirus B` + 0.665 * `Enterovirus C` + 0.708 * `Enterovirus D` + 0.538 * `Rhinovirus A` + 0.571 * `Rhinovirus B`)/5) %>% 
  mutate(`EVB adjusted` = `Enterovirus B` - (0.789 * `Enterovirus A` + 0.773 * `Enterovirus C` + 0.732 * `Enterovirus D` + 0.592 * `Rhinovirus A` + 0.573 * `Rhinovirus B`)/5) %>% 
  mutate(`EVC adjusted` = `Enterovirus C` - (0.773 * `Enterovirus B` + 0.665 * `Enterovirus A` + 0.755 * `Enterovirus D` + 0.509 * `Rhinovirus A` + 0.626 * `Rhinovirus B`)/5) %>% 
  mutate(`EVD adjusted` = `Enterovirus D` - (0.732 * `Enterovirus B` + 0.755 * `Enterovirus C` + 0.708 * `Enterovirus A` + 0.717 * `Rhinovirus A` + 0.722 * `Rhinovirus B`)/5) %>% 
  mutate(`RVA adjusted` = `Rhinovirus A` - (0.592 * `Enterovirus B` + 0.509 * `Enterovirus C` + 0.717 * `Enterovirus D` + 0.538 * `Enterovirus A` + 0.656 * `Rhinovirus B`)/5) %>% 
  mutate(`RVB adjusted` = `Rhinovirus B` - (0.573 * `Enterovirus B` + 0.626 * `Enterovirus C` + 0.722 * `Enterovirus D` + 0.656 * `Rhinovirus A` + 0.571 * `Enterovirus A`)/5)

ggpairs(data = HSV1_HSV2_deconv,
        columns = c(9:14),
        #mapping = aes(alpha = .1)
        lower = list(continuous = pairs_dots)) +
  theme_light() + 
  theme(text = element_text(size = 8))
