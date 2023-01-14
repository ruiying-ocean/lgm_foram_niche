# grouping EPILOG data
library(data.table)
library(tidyverse)

## read symbiosis table
source("code/data_clean/read_symbiosis_table.R")

## read and select data (remove unidentified species as well)
glamap <- read_csv("data/lgm_foram_count_data/LGM_GLAMAP.csv")
glamap <- glamap %>% pivot_longer(cols=c(`G. aequilateralis [%]`:`N. pachyderma d [%]`), names_to = "Species", values_to="Relative_Abundance")
glamap <- glamap %>%  mutate_at("Species", str_replace, " \\[%\\]", "")

glamap <- glamap %>% mutate(Species=recode(Species, 
                                           "N. pachyderma d" =  "N. pachyderma",
                                           "N. pachyderma s" = "N. incompta",
                                           "G. ruber w" = "G. ruber",
                                           "G. ruber p" = "G. ruber",
                                           "G. ruber hsp" = "G. ruber",
                                           "G. quinqueloba d"= "T. quinqueloba",
                                           "G. quinqueloba" = "T. quinqueloba",
                                           "G. quinqueloba s" =  "T. quinqueloba",
                                           "G. truncatulinoides s" = "G. truncatulinoides",
                                           "G. truncatulinoides d" = "G. truncatulinoides",                                           
                                           "P/D int...44"="N. dutertrei",
                                           "P/D int...26"="N. dutertrei",
                                           "G. trilobus sac" = "G. sacculifer",
                                           "G. trilobus tril" = "G. clavaticamerata",
                                           "G. mentum" = "O. NA",
                                           "D. grahami" = "O.NA"
))


find_sp(glamap)

glamap_merged <- merge(glamap, symbiosis_tbl, by.x="Species", by.y = "short_name") %>% select(!c(Species, Species.y))

# aggregate functional group abundance and divided by 100
glamap_merged <- glamap_merged %>% group_by(Campaign, Event, Latitude, Longitude, `Date/Time`, `Depth [m]`, `Elevation [m]`,  Symbiosis, Spinose) %>% 
    summarise_all(.funs = sum, na.rm=T) %>% ungroup() %>% filter(Symbiosis != "Undetermined" & Spinose != 'Undetermined')

glamap_merged <- glamap_merged %>% group_by(Campaign, Event, Latitude, Longitude, `Date/Time`, `Elevation [m]`,  Symbiosis, Spinose) %>% 
    summarise_all(.funs = mean, na.rm=T) %>% ungroup()

glamap_merged <- glamap_merged %>% mutate(Relative_Abundance = Relative_Abundance/100)

# export
write_csv(glamap_merged, "data/lgm_foram_count_data/LGM_GLAMAP_groupped.csv")
