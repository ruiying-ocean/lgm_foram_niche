# grouping EPILOG data
library(data.table)
library(tidyverse)

## read symbiosis table
source("code/data_clean/read_symbiosis_table.R")

## read and calculate relative abundance
climap <- read_csv("data/lgm_foram_count_data/LGM_CLIMAP.csv")
climap <- climap %>% mutate(total = rowSums(across(`O. universa [#]`:`Foram plankt oth [#]`))) %>%
  mutate(across(`O. universa [#]`:`Foram plankt oth [#]`, ~ ./ total)) %>% select(-total)

# clean species name
climap <- climap %>% pivot_longer(cols=c(`O. universa [#]`:`Foram plankt oth [#]`), names_to = "Species", values_to="Relative_Abundance")
climap <- climap %>%  mutate_at("Species", str_replace, " \\[#\\]", "")

climap <- climap %>% mutate(Species=recode(Species, 
                                           "N. pachyderma d" =  "N. pachyderma",
                                           "G. pachyderma" = "N. pachyderma",
                                           "G. ruber (sum of G. ruber pink and G. r...)" ="G. ruber",
                                           "G. anfracta" = "D. anfracta",
                                           "N. pachyderma s" = "N. incompta",
                                           "G. ruber w" = "G. ruber",
                                           "G. ruber p" = "G. ruber",
                                           "G. ruber hsp" = "G. ruber",
                                           "G. quinqueloba d"= "T. quinqueloba",
                                           "G. quinqueloba" = "T. quinqueloba",
                                           "G. quinqueloba s" =  "T. quinqueloba",
                                           "G. truncatulinoides s" = "G. truncatulinoides",
                                           "G. truncatulinoides d" = "G. truncatulinoides",                                           
                                           "G. trilobus sac" = "G. sacculifer",
                                           "G. sacculifer sac" = "G. sacculifer",
                                           "G. sacculifer wo sac" = "G. sacculifer",
                                          "G. humilis" = "T. humilis",
                                          "G. hexagona" = "G. hexagonus",
                                           "G. sacculifer (sum of G. sacculifer no sac a...)" = "G. sacculifer",
                                           "G. tumida flexuosa" = "G. tumida",
                                           "G. pumilio"= "B. pumilio",
                                           "G. iota" = "T. iota",
                                           "G. dutertrei" ="N. dutertrei",
                                           "G. trilobus tril" = "G. clavaticamerata",
                                           "G. menardii (sum of G. menardii, G. tumida...)"="G. menardii",
                                           "G. bradyi" = "G. uvula",
                                           "Foram plankt oth" = "O. NA"
))


find_sp(climap)

climap_merged <- merge(climap, symbiosis_tbl, by.x="Species", by.y = "short_name") %>% select(!c(Species, Species.y))

# aggregate functional group abundance and divided by 100
climap_merged <- climap_merged %>% group_by(Latitude, Longitude, Event, `Depth [m]`, `Elevation [m]`,  Symbiosis, Spinose) %>% 
    summarise_all(.funs = sum, na.rm=T) %>% ungroup() %>% filter(Symbiosis != "Undetermined" & Spinose != 'Undetermined')

climap_merged <- climap_merged %>% group_by(Latitude, Longitude, Event, `Elevation [m]`,  Symbiosis, Spinose) %>% 
    summarise_all(.funs = mean, na.rm=T) %>% ungroup()

write_csv(climap_merged, "data/lgm_foram_count_data/LGM_CLIMAP_groupped.csv")
