## This script contains ForamEcoGENIE output
## Contact: rui.ying@bristol.ac.uk

source('code/lib.R')

## Model output
genie_fg_raw <- load_models("model/model_drived/")
## only get columns in the pattern 'xx_c'

## carbon biomass
## species here is actually ecogroup
genie_fg_raw <- genie_fg_raw %>% select(-c("sn")) %>% 
  pivot_longer(cols=bn:ss, names_to = "species", values_to="biomass") %>%
  convert_to_abundance() %>%
  mutate(species=recode(species,
                          "bn"="Symbiont-barren Non-Spinose",
                          "bs"="Symbiont-barren Spinose",
                          "ss"="Symbiont-obligate Spinose"))

quantlvl <- seq(0.9, 0.99, 0.01)

genie_fg_smooth <- loop_smooth(genie_fg_raw, i = species, j = age, x=sst, y=abundance_michaels, quant_level=quantlvl)

## save in Rdata
save(genie_fg_smooth, file = "data/genie_fg_smooth.Rdata")
save(genie_fg_raw, file = "data/genie_fg_raw.Rdata")

### -----------------------
### basic statistics

### report total data points
## obs_fg_a_raw %>% group_by(species, age) %>% summarise(n())
## obs_sp_raw %>% group_by(species, age) %>% summarise(n())


## ------------------------
## trait does not explain the difference

# sp_analaysis <- thermal_opt(obs_sp_smooth) %>% group_by(age, species) %>% summarise(T_opt =model_x) %>%
#   pivot_wider(id_cols = "species",names_from = "age",values_from = "T_opt") %>% mutate(diff=LGM-PI)
# sp_analaysis %>% write_csv("data/model_drived/Topt_sp_lgm.csv")
# trait_info <- read.csv("~/Science/lgm_foram_census/fg/foram_sp_db.csv") %>%
#   mutate(sp = map_vec(Name, species_abbrev)) %>% select(sp, Symbiosis, Spinose)
# 
# sp_analaysis <- merge(sp_analaysis,trait_info, by.x="species",by.y="sp")
# aov(diff ~ Symbiosis * Spinose, data=sp_analaysis) %>% summary()
