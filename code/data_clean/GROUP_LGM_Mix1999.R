# grouping EPILOG data
library(data.table)
library(tidyverse)

## read symbiosis table
source("code/data_clean/read_symbiosis_table.R")

## read and select data (remove unidentified species as well)
mix1999 <- read_csv("data/lgm_foram_count_data/LGM_MIX1999.csv")
mix1999 <- mix1999 %>% pivot_longer(cols=c(`O. universa [%]`:`G. glutinata [%]`), names_to = "Species", values_to="Relative_Abundance")
mix1999 <- mix1999 %>%  mutate_at("Species", str_replace, " \\[%\\]", "")

mix1999 <- mix1999 %>% mutate(Species=recode(Species, 
                                             "N. pachyderma d" =  "N. pachyderma",
                                             "N. pachyderma s" = "N. incompta",     
                                             "G. quinqueloba"="T. quinqueloba",
                                             "G. hexagona" = "G. hexagonus",                                             
                                           "G. truncatulinoides s" = "G. truncatulinoides",
                                           "G. truncatulinoides d" = "G. truncatulinoides"))

## check every species is in the list
find_sp(mix1999)

# aggregate functional group's abundance
mix1999_merged <- merge(mix1999, symbiosis_tbl, by.x="Species", by.y = "short_name") %>% select(!c(Species, Species.y))

mix1999_merged <- mix1999_merged %>% group_by(Event, Latitude, Longitude, `Elevation [m]`, `Depth [m]`,`SST (1-12) [°C]`, Symbiosis, Spinose) %>% 
    summarise_all(.funs = sum, na.rm=T) %>% ungroup() %>% filter(Symbiosis != "Undetermined" & Spinose != 'Undetermined')

mix1999_merged <- mix1999_merged %>% group_by(Event, Latitude, Longitude, `Elevation [m]`, `SST (1-12) [°C]`, Symbiosis, Spinose) %>% 
  summarise_all(.funs = mean, na.rm=T) %>% ungroup()

mix1999_merged <- mix1999_merged %>% mutate(Relative_Abundance = Relative_Abundance/100)

write_csv(mix1999_merged, "data/lgm_foram_count_data/LGM_MIX1999_groupped.csv")
