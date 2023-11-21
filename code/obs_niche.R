## This script contains observational thermal performance data
## Contact: rui.ying@bristol.ac.uk

source("code/lib.R") ## load functions

## Read raw data
## OBSERVATION DATA (absolute and functional group format)
lgm_fg_a <- read_csv("https://raw.githubusercontent.com/ruiying-ocean/lgm_foram_census/main/tidy/lgm_fg_a_wsst.csv") %>%
  mutate(age = "LGM") %>%
  select(c("Latitude", "Longitude", "age", "SST", "symbiont-barren non-spinose", "symbiont-barren spinose", "symbiont-obligate spinose"))
lgm_fg_a <- lgm_fg_a %>% pivot_longer(cols = `symbiont-barren non-spinose`:`symbiont-obligate spinose`, names_to = "species", values_to = "abundance")

pi_fg_a <- read_csv("https://raw.githubusercontent.com/ruiying-ocean/lgm_foram_census/main/tidy/forcens_fg_a_wsst.csv") %>%
  mutate(age = "PI") %>%
  select(c("Latitude", "Longitude", "age", "SST", "symbiont-barren non-spinose", "symbiont-barren spinose", "symbiont-obligate spinose"))
pi_fg_a <- pi_fg_a %>% pivot_longer(cols = `symbiont-barren non-spinose`:`symbiont-obligate spinose`, names_to = "species", values_to = "abundance")

lgm_sp <- read_csv("https://raw.githubusercontent.com/ruiying-ocean/lgm_foram_census/main/tidy/lgm_sp_a_wsst.csv") %>% mutate(age = "LGM")
lgm_sp <- lgm_sp %>% pivot_longer(cols = `G. bulloides`:`G. elongatus`, names_to = "species", values_to = "Abundance")

pi_sp <- read_csv("https://raw.githubusercontent.com/ruiying-ocean/lgm_foram_census/main/tidy/forcens_sp_a_wsst.csv") %>% mutate(age = "PI")
pi_sp <- pi_sp %>% pivot_longer(cols = `D. anfracta`:`H. digitata`, names_to = "species", values_to = "Abundance")

## combine ecogroup data frames and
## exclude the symbiont-facultative non-spinose ecogroup
obs_fg_a_raw <- rbind(pi_fg_a, lgm_fg_a) %>% rename(sst = SST)
obs_fg_a_raw <- obs_fg_a_raw %>% filter(species != "symbiont-facultative non-spinose")

## only get relevant columns from the observation data (species-level)
subset_columns <- c("SST", "species", "Abundance", "age")
obs_sp_raw <- rbind(pi_sp[subset_columns], lgm_sp[subset_columns])

## exclude species with too low abundance
lgm_sp_list <- lgm_sp %>%
  pull(species) %>%
  unique() ## 44 spcies
pi_sp_list <- pi_sp %>%
  pull(species) %>%
  unique() ## 41 species

sp_list <- intersect(lgm_sp_list, pi_sp_list) ## 37 species

## exclude species with too low abundance
sp_list <- sp_list[!sp_list %in% c(
  "G. uvula", "H. digitata", "B. pumilio"
)]

obs_sp_raw <- obs_sp_raw %>% filter(species %in% sp_list)

## Smooth the data
quantlvl <- seq(0.9, 0.99, 0.01)
obs_fg_a_smooth <- loop_smooth(obs_fg_a_raw, i = species, j = age, x = sst, y = abundance, quant_level = quantlvl)
obs_sp_smooth <- loop_smooth(obs_sp_raw, i = species, j = age, x = SST, y = Abundance, quant_level = quantlvl)

## export the data in Rdata
save(obs_fg_a_smooth, obs_sp_smooth, file = "data/obs_smooth.Rdata")
save(obs_fg_a_raw, obs_sp_raw, file = "data/obs_raw.Rdata")
