# grouping EPILOG data
library(data.table)
library(tidyverse)

## read symbiosis table
source("code/data_clean/read_symbiosis_table.R")

## read and select data (remove unidentified species as well)
glamap <- read_csv("data/modern_data/GLAMAP_modern.csv")
glamap <- glamap %>% pivot_longer(cols=c(`G. bulloides [%]`:`G. menardii [%]`), names_to = "Species", values_to="Relative_Abundance")
glamap <- glamap %>%  mutate_at("Species", str_replace, " \\[%\\]", "")

glamap <- glamap %>% mutate(Species=recode(Species, 
                                           "N. pachyderma d" =  "N. pachyderma",
                                           "N. pachyderma s" = "N. incompta",
                                           "G. ruber w" = "G. ruber",
                                           "G. ruber p" = "G. ruber",
                                           "G. ruber hsp" = "G. ruber",
                                           "G. sacculifer wo sac" = "G. sacculifer" ,
                                           "G. sacculifer sac" = "G. sacculifer" ,
                                           "G. quinqueloba d"= "T. quinqueloba",
                                           "G. quinqueloba" = "T. quinqueloba",
                                           "G. bradyi" = "G. uvula",
                                           "G. quinqueloba s" =  "T. quinqueloba",
                                           "G. truncatulinoides s" = "G. truncatulinoides",
                                           "G. truncatulinoides d" = "G. truncatulinoides",                                           
                                           "G. trilobus sac" = "G. sacculifer",
                                           "G. trilobus tril" = "G. clavaticamerata",
))

find_sp(glamap)

glamap_merged <- merge(glamap, symbiosis_tbl, by.x="Species", by.y = "short_name") %>% select(!c(Species, Species.y))

# aggregate functional group abundance and divided by 100
glamap_merged <- glamap_merged %>% group_by(Event, Latitude, Longitude, `Depth [m]`, `Elevation [m]`,  Symbiosis, Spinose) %>% 
    summarise_all(.funs = sum, na.rm=T) %>% ungroup() %>% filter(Symbiosis != "Undetermined" & Spinose != 'Undetermined')

glamap_merged <- glamap_merged %>% group_by(Event, Latitude, Longitude, `Elevation [m]`,  Symbiosis, Spinose) %>% 
    summarise_all(.funs = mean, na.rm=T) %>% ungroup()

glamap_merged <- glamap_merged %>% mutate(Relative_Abundance = Relative_Abundance/100)

write_csv(glamap_merged, "data/modern_data/GLAMAP_modern_groupped.csv")

