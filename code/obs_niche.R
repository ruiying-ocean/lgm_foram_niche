## This script contains observational thermal performance data
## Contact: rui.ying@bristol.ac.uk

source("code/lib.R") ## load functions

## Read raw data
## OBSERVATION DATA (absolute and functional group format)
# dir <- "https://raw.githubusercontent.com/ruiying-ocean/lgm_foram_census/main/tidy"
dir <- "~/Science/lgm_foram_census/tidy"
lgm_fg_r <- read_csv(paste(dir, "lgm_fg_r_wsst.csv", sep = "/"))

lgm_fg_r <- lgm_fg_r %>%
  mutate(age = "LGM") %>%
  dplyr::filter(Data_Source == "margo") %>%
  select(c("Latitude", "Longitude", "age", "SST", "symbiont-barren non-spinose", "symbiont-barren spinose", "symbiont-obligate spinose"))

lgm_fg_r <- lgm_fg_r %>% pivot_longer(cols = `symbiont-barren non-spinose`:`symbiont-obligate spinose`, names_to = "species", values_to = "abundance")

pi_fg_r <- read_csv(paste(dir, "forcens_fg_r_wsst.csv", sep = "/")) %>%
  mutate(age = "PI") %>%
  select(c("Latitude", "Longitude", "age", "SST", "symbiont-barren non-spinose", "symbiont-barren spinose", "symbiont-obligate spinose"))

pi_fg_r <- pi_fg_r %>% pivot_longer(cols = `symbiont-barren non-spinose`:`symbiont-obligate spinose`, names_to = "species", values_to = "abundance")

lgm_sp_r <- read_csv(paste(dir, "lgm_sp_r_wsst.csv", sep = "/")) %>%
  mutate(age = "LGM") %>%
  dplyr::filter(Data_Source == "margo")

lgm_sp_r <- lgm_sp_r %>% pivot_longer(cols = `O. universa`:`G. siphonifera`, names_to = "species", values_to = "Abundance")

pi_sp_r <-read_csv(paste(dir, "forcens_sp_r_wsst.csv", sep = "/")) %>%
  mutate(age = "PI")

pi_sp_r <- pi_sp_r %>% pivot_longer(cols = `D. anfracta`:`H. digitata`, names_to = "species", values_to = "Abundance")

## combine ecogroup data frames and
## exclude the symbiont-facultative non-spinose ecogroup
obs_fg_r_raw <- rbind(pi_fg_r, lgm_fg_r) %>% rename(sst = SST)
obs_fg_r_raw <- obs_fg_r_raw %>% dplyr::filter(species != "symbiont-facultative non-spinose")

## only get relevant columns from the observation data (species-level)
subset_columns <- c("SST", "species", "Abundance", "age")
obs_sp_r_raw <- rbind(pi_sp_r[subset_columns], lgm_sp_r[subset_columns])

lgm_sp_list <- lgm_sp_r %>%
    dplyr::filter(Abundance > 0) %>%
    pull(species) %>%
    unique()

pi_sp_list <- pi_sp_r %>%
    dplyr::filter(Abundance > 0) %>%
    pull(species) %>%
    unique()

sp_list <- intersect(lgm_sp_list, pi_sp_list)

obs_sp_r_raw <- obs_sp_r_raw %>% dplyr::filter(species %in% sp_list) %>% drop_na(Abundance)

## Smooth the data
quantlvl <- seq(0.9, 0.99, 0.01)
obs_fg_r_smooth <- loop_smooth(obs_fg_r_raw, i = species, j = age, x = sst, y = abundance, quant_level = quantlvl)
obs_sp_r_smooth <- loop_smooth(obs_sp_r_raw, i = species, j = age, x = SST, y = Abundance, quant_level = quantlvl)

## export the data in Rdata
save(obs_fg_r_smooth, obs_sp_r_smooth, file = "data/obs_smooth.Rdata")
save(obs_fg_r_raw, obs_sp_r_raw, file = "data/obs_raw.Rdata")

## export the statistic data in csv
obs_sp_r_smooth %>% thermal_opt(long_format=F) %>%
    mutate(Topt_mean_diff=PI_Topt_mean-LGM_Topt_mean) %>%
    write_csv("data/Topt_sp_lgm.csv")

obs_fg_r_smooth %>% thermal_opt(long_format=F) %>%
    mutate(Topt_mean_diff=PI_Topt_mean-LGM_Topt_mean) %>%
    write_csv("data/Topt_fg_lgm.csv")
