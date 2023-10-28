## This script contains ForamEcoGENIE output
## Contact: rui.ying@bristol.ac.uk


## Model output
genie_fg_raw <- load_models("data/model_drived/")
## only get columns in the pattern 'xx_c'

## carbon biomass
## species here is actually ecogroup
genie_fg_c <- genie_fg_raw %>%
  select(-contains("p"), -contains("fe")) %>%
  pivot_longer(cols=bn_c:ss_c, names_to = "species", values_to="biomass") %>%
  convert_to_abundance() %>%
  mutate(species=recode(species,
                          "bn"="Symbiont-barren Non-Spinose",
                          "bs"="Symbiont-barren Spinose",
                          "ss"="Symbiont-obligate Spinose"))

## convert c(0.9, 0.95) to c("model_y_0.9", "model_y_0.95")



genie_fg_smooth <- loop_smooth(genie_fg_c, i = species, j = age, x=sst, y=abundance_michaels, quant_level=quantlvl)
genie_fg_smooth_chl <- loop_smooth(genie_fg_c, i = species, j = age, x=chl_total, y=abundance_michaels, quant_level=quantlvl)
genie_fg_smooth_light <- loop_smooth(genie_fg_c, i = species, j = age, x=light, y=abundance_michaels, quant_level=quantlvl)

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
