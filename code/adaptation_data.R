library(tidyverse)

## 1. abundance and temperature data

## -------- Taxonomy ------------
## G. sacculifer with or without sac-like final chamber are merged
## Pachederma/Dutertrei integrade is not considered
## G. menardii is considered as G. cultrata
## G. Tumida, G. menardii flexuosa are merged as G. flexuosa
## G. ruber (w) and G. ruber (p) are separated to G. ruber albus (white in latin) and G. ruber ruber (reddish in latin)

select_species <- c("Orbulina universa", "Globigerina bulloides", "Neogloboquadrina pachyderma",
  "Neogloboquadrina dutertrei", "Neogloboquadrina incompta", "Globorotalia inflata", "Globigerinita glutinata",
  "Globigerinoides ruber w", "Globigerinoides ruber p", "Globorotalia menardii" ,
  "Turborotalita quinqueloba","Trilobatus sacculifer"
  )

## T. quinqueloba upper-thermocline or mixed-layer habitat
## T. sacculifer  Open ocean mixed-layer tropical/subtropical
## G. menardii Open ocean thermocline
## G. glutinata  shallow, mixed-layer habitat
## G. inflata Open ocean thermocline.

lgm_climap <- read_csv("data/lgm_foram/raw/LGM_CLIMAP_wsst.csv")
lgm_atlantic <- read_csv("data/lgm_foram/raw/LGM_MARGO/LGM_MARGO_SouthAtlantic_count_wsst.csv")
lgm_pacific <- read_csv("data/lgm_foram/raw/LGM_MARGO/LGM_MARGO_pacific_count_wsst.csv")

names(lgm_climap) <- gsub(" [#]", "", names(lgm_climap), fixed=T)
names(lgm_climap) <- gsub(" [m]", "", names(lgm_climap), fixed=T)

#`G. truncatulinoides` is strangely high          
lgm_climap <- lgm_climap %>%
  select(SST,SSS,`O. universa`, `G. bulloides`,`G. inflata`,`G. quinqueloba`,
         `N. dutertrei`, `G. ruber w`, `G. ruber p`,
         `N. pachyderma d`, `N. pachyderma s`, `G. glutinata`,
         `G. menardii`,
         `G. sacculifer total`
         ) %>%
  rename(`N. incompta` = `N. pachyderma d`,
         `N. pachyderma` = `N. pachyderma s`,
         `T. quinqueloba`=`G. quinqueloba`,
         `T. sacculifer` = `G. sacculifer total`)
  
lgm_atlantic <- lgm_atlantic %>% rename(`Neogloboquadrina pachyderma` = `Neogloboquadrina pachyderma L`,
                        `Neogloboquadrina incompta`= `Neogloboquadrina pachyderma R`,
                        `Globigerinita glutinata` = `Globigerinita glutinata`,
                        `Globigerinoides ruber p` = `Globigerinoides ruber (pink)`,
                        `Globigerinoides ruber w` = `Globigerinoides ruber (white)`,
                        `Trilobatus sacculifer`  = `Globigerinoides sacc total`) %>%
  select(SST,SSS, all_of(select_species))

lgm_pacific <- lgm_pacific %>%
  rename(`Neogloboquadrina pachyderma` = `Neogloboquadrina pachyderma L`,
         `Neogloboquadrina incompta`= `Neogloboquadrina pachyderma R`,
         `Globigerinita glutinata` = `Globigerinita glutinata`,) %>%
  select(SST,SSS, any_of(select_species))

### G ruber pink and white!
species_abbrev <- function(full_name, split_str = " ", sep_string=". "){
  splited_string <- str_split(full_name, split_str)[[1]]
  
  genus_name <- splited_string[1]
  sp_name <- splited_string[2]
  
  genus_abbrev <- str_sub(genus_name, 1, 1)
  combine_name <- paste(genus_abbrev, sp_name, sep = sep_string)
  
  if (length(splited_string) > 2){
    others <- str_split(full_name, split_str)[[1]][3]
    all_name <- paste0(combine_name, " ", others)
    return(all_name)
  } else{
    return(combine_name)
  }
}

lgm_pacific <- lgm_pacific %>% pivot_longer(cols =  !c(SST,SSS), names_to = "Species", values_to = "Count") %>%
  mutate(Species=sapply(Species, species_abbrev))

lgm_atlantic <- lgm_atlantic %>% pivot_longer(cols = !c(SST,SSS), names_to = "Species", values_to = "Count") %>%
  mutate(Species=sapply(Species, species_abbrev))

lgm_climap <- lgm_climap %>% pivot_longer(cols =  !c(SST,SSS), names_to = "Species", values_to = "Count") %>%
  mutate(Species=sapply(Species, species_abbrev))

lgm_niche <- rbind(lgm_pacific, lgm_atlantic, lgm_climap)
lgm_niche$age = "lgm"

### ForCenS
pi <- read_csv("data/modern_foram/raw/forcens_raw_count_wsst.csv")
pi_niche <-pi %>% select(c("Orbulina_universa", "Globigerina_bulloides", "Neogloboquadrina_pachyderma",
                    "Neogloboquadrina_dutertrei", "Neogloboquadrina_incompta", "Globoconella_inflata",
                    "Globigerinita_glutinata", "Turborotalita_quinqueloba",
                    "Trilobatus_sacculifer", "Globorotalia_menardii",
                    "Globigerinoides_white", "Globigerinoides_ruber"), SST, SSS) %>%
  rename(Globigerinoides_ruber_w=Globigerinoides_white,
         Globigerinoides_ruber_p=Globigerinoides_ruber)
pi_niche <- pi_niche %>%
  pivot_longer(cols = !c(SST,SSS), names_to = "Species", values_to = "Count") %>%
  mutate(Species=sapply(Species, species_abbrev, split_str="_"))

pi_niche$age = "pi"
#source("code/lombard_2009_model.R")

## Modern sediment trap data
modern_niche <- read_csv("data/Jonkers_2019/shell_flux_data_wsst.csv")
names(modern_niche)[20:60] <- sapply(names(modern_niche)[20:60], species_abbrev,  split_str="_")
modern_niche$SSS <- NA
modern_niche <- modern_niche %>% rename(`G. ruber w` = `G. white`, 
                                        `G. ruber p` = `G. ruber`,
                                        #`G. quinqueloba`=`T. quinqueloba`,
                                        `G. menardii` = `G. merdii`,
                                        `G. glutinata`=`G. glutita`)
modern_niche <- modern_niche %>%  select(SST, SSS, `O. universa`, `G. bulloides`,`G. inflata`,`T. quinqueloba`,
       `N. dutertrei`, `G. ruber w`, `G. ruber p`,
       `N. pachyderma`, `N. incompta`, `G. glutinata`,
       `G. menardii`, `T. sacculifer`)

modern_niche <- modern_niche %>% pivot_longer(cols=`O. universa`:`T. sacculifer`, names_to = "Species", values_to = "Flux")

# flux is in ind m-2 day-1
modern_niche <- modern_niche %>% rename(Count=Flux) %>% drop_na(SST)
modern_niche$age = "modern"
modern_niche %>% mutate(Count=Count/1000)
niche_data <- rbind(pi_niche,lgm_niche) %>% drop_na(SST)
niche_data <- niche_data %>% rename(sst=SST, sss=SSS, count=Count, species=Species)

remove(lgm_atlantic)
remove(lgm_climap)
remove(lgm_niche)
remove(lgm_pacific)
remove(modern_niche)
remove(pi_niche)
remove(pi)

## 2. biomass and temperature model output

model_modern <- read_csv("data/modern_4pca.csv") %>% select(-1) %>% mutate(age  = "modern")
model_lgm <- read_csv("data/lgm_4pca.csv") %>% select(-1)  %>% mutate(age  = "LGM")
niche_genie <- rbind(model_modern, model_lgm) %>% pivot_longer(cols=bn:ss, names_to = "species", values_to="biomass")
remove(model_modern)
remove(model_lgm)
