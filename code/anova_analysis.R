## This script is to analyse the species optimal temperature change
## and its behind driver, here as the symbiont and spine trait

## Contact: rui.ying@bristol.ac.uk

## load visualisation packages
source("code/lib.R")
library(ggpubr)

## read presaved Rdata
load("data/obs_smooth.Rdata")
load("data/genie_fg_smooth.Rdata")

##  as the Fig. 2
## delta habitat temperature
delta_habtemp <- read_csv("data/obs_sp_hab.csv")
obs_sp_Topt <- obs_sp_Topt %>% left_join(delta_habtemp)

## merge with species trait data table
foram_sp_db <- read_csv("https://raw.githubusercontent.com/ruiying-ocean/lgm_foram_census/main/fg/foram_taxonomy.csv")

trait_info <-foram_sp_db %>%
   mutate(sp = map_vec(`Species name`, species_abbrev)) %>% select(sp, Symbiosis, Spine)

## check all species are included
common_sp <- intersect(obs_sp_Topt$species, trait_info$sp)
if (length(common_sp) == length(obs_sp_Topt$species)) {
  print("all species included")
}

obs_sp_Topt <- obs_sp_Topt %>% rename(sp=species) %>% left_join(trait_info, by = "sp") %>% distinct()

## use species trait to explain the species optimal temperature change
mod <- aov(Topt_mean_diff ~ Spine + Symbiosis, data = obs_sp_Topt)
summary(mod)

## use species optimal temperature to explain the species optimal temperature change
ext_data_fig6b <- plot_lm(data=obs_sp_Topt, x="PI_Topt_mean", y="Topt_mean_diff", c(9,6,13,5))
ext_data_fig6b <- ext_data_fig6b + theme_publication(base_size = 14) + 
  labs(x = "PI habitat mean temperature (°C)", y = "Δ Thermal optimum (°C)")

## use species habitat temperature change to explain the species optimal temperature change
ext_data_fig6a <- plot_lm(data=obs_sp_Topt, x="delta_habtemp", y="Topt_mean_diff", c(3,-4, 4.2,-5))
ext_data_fig6a <- ext_data_fig6a + theme_publication(base_size = 14) + 
  labs(x = "Habitat temperature change (°C)", y = "Δ Thermal optimum (°C)")

library(patchwork)
ext_data_fig6 <- ext_data_fig6a + ext_data_fig6b + plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a")
ext_data_fig6 %>% ggsave(file = "output/ext_data_fig6.png", dpi = 300, width = 10, height = 4)
