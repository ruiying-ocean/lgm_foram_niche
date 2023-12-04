## This script contains ForamEcoGENIE output
## Contact: rui.ying@bristol.ac.uk

source('code/lib.R')

## Model output
genie_fg_raw <- load_models("model/model_drived/")

## carbon biomass
## species here is actually ecogroup
genie_fg_raw <- genie_fg_raw %>% select(-c("sn")) %>% 
  pivot_longer(cols=bn:ss, names_to = "species", values_to="abundance") %>%
  mutate(species=recode(species,
                          "bn"="symbiont-barren non-spinose",
                          "bs"="symbiont-barren spinose",
                          "ss"="symbiont-obligate spinose"))

quantlvl <- seq(0.9, 0.99, 0.01)

genie_fg_smooth <- loop_smooth(genie_fg_raw, i = species, j = age, x=sst, y=abundance, quant_level=quantlvl)

genie_fg_Topt <- genie_fg_smooth %>%  thermal_opt(long_format=F) %>% mutate_if(is.numeric, round, 1)

## save in Rdata
save(genie_fg_smooth, genie_fg_Topt, file = "data/genie_fg_smooth.Rdata")
save(genie_fg_raw, file = "data/genie_fg_raw.Rdata")

