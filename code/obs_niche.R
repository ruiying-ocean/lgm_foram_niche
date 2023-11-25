## This script contains observational thermal performance data
## Contact: rui.ying@bristol.ac.uk

source("code/lib.R") ## load functions

## Read raw data
## OBSERVATION DATA (absolute and functional group format)
lgm_fg_r <- read_csv("https://raw.githubusercontent.com/ruiying-ocean/lgm_foram_census/main/tidy/lgm_fg_r_wsst.csv")

lgm_fg_r <- lgm_fg_r %>%
  mutate(age = "LGM") %>%
  dplyr::filter(Data_Source == "margo") %>%
  select(c("Latitude", "Longitude", "age", "SST", "symbiont-barren non-spinose", "symbiont-barren spinose", "symbiont-obligate spinose"))

lgm_fg_r <- lgm_fg_r %>% pivot_longer(cols = `symbiont-barren non-spinose`:`symbiont-obligate spinose`, names_to = "species", values_to = "abundance")

pi_fg_r <- read_csv("https://raw.githubusercontent.com/ruiying-ocean/lgm_foram_census/main/tidy/forcens_fg_r_wsst.csv") %>%
  mutate(age = "PI") %>%
  select(c("Latitude", "Longitude", "age", "SST", "symbiont-barren non-spinose", "symbiont-barren spinose", "symbiont-obligate spinose"))

pi_fg_r <- pi_fg_r %>% pivot_longer(cols = `symbiont-barren non-spinose`:`symbiont-obligate spinose`, names_to = "species", values_to = "abundance")

lgm_sp_r <- read_csv("https://raw.githubusercontent.com/ruiying-ocean/lgm_foram_census/main/tidy/lgm_sp_r_wsst.csv") %>%
  mutate(age = "LGM") %>%
  dplyr::filter(Data_Source == "margo")
lgm_sp_r <- lgm_sp_r %>% pivot_longer(cols = `O. universa`:`G. siphonifera`, names_to = "species", values_to = "Abundance")

pi_sp_r <- read_csv("https://raw.githubusercontent.com/ruiying-ocean/lgm_foram_census/main/tidy/forcens_sp_r_wsst.csv") %>% mutate(age = "PI")
pi_sp_r <- pi_sp_r %>% pivot_longer(cols = `D. anfracta`:`H. digitata`, names_to = "species", values_to = "Abundance")

## combine ecogroup data frames and
## exclude the symbiont-facultative non-spinose ecogroup
obs_fg_r_raw <- rbind(pi_fg_r, lgm_fg_r) %>% rename(sst = SST)
obs_fg_r_raw <- obs_fg_r_raw %>% dplyr::filter(species != "symbiont-facultative non-spinose")

## only get relevant columns from the observation data (species-level)
subset_columns <- c("SST", "species", "Abundance", "age")
obs_sp_r_raw <- rbind(pi_sp_r[subset_columns], lgm_sp_r[subset_columns])
## exclude species with too low abundance
lgm_sp_list <- lgm_sp_r %>%
  pull(species) %>%
  unique()
pi_sp_list <- pi_sp_r %>%
  pull(species) %>%
  unique()

sp_list <- intersect(lgm_sp_list, pi_sp_list)

## exclude species with low abundance
sp_list <- sp_list[!sp_list %in% c(
  "G. uvula", "H. digitata", "B. pumilio",
  "T. iota","H. pelagica", "D. anfracta",
  "G. adamsi", "G. eastropacia"
)]

obs_sp_r_raw <- obs_sp_r_raw %>% dplyr::filter(species %in% sp_list) %>% drop_na(Abundance)

## Smooth the data
quantlvl <- seq(0.9, 0.99, 0.01)
obs_fg_r_smooth <- loop_smooth(obs_fg_r_raw, i = species, j = age, x = sst, y = abundance, quant_level = quantlvl)
obs_sp_r_smooth <- loop_smooth(obs_sp_r_raw, i = species, j = age, x = SST, y = Abundance, quant_level = quantlvl)

## export the data in Rdata
save(obs_fg_r_smooth, obs_sp_r_smooth, file = "data/obs_smooth.Rdata")
save(obs_fg_r_raw, obs_sp_r_raw, file = "data/obs_raw.Rdata")
