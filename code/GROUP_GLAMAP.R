# grouping EPILOG data
library(data.table)
library(tidyverse)

## read symbiosis table
symbiosis_tbl <- fread('~/ForamData/data/Symbiosis_table.csv')
symbiosis_tbl[, "Name" := lapply(.SD, function(x) gsub(" ", "_", x)),
              .SDcol="Name"] #replace white space by underscore(_)
symbiosis_tbl <- symbiosis_tbl[, -("Remark")] #delete comment column
symbiosis_tbl <- rbindlist(list(symbiosis_tbl, data.table(Name="Others",
                                                          Spinose="Undetermined",
                                                          Symbiosis="Undetermined")),
                           use.names = TRUE)
setnames(symbiosis_tbl, "Name", "Species")
symbiosis_tbl$Species <- gsub("_", " ",  symbiosis_tbl$Species)

## Add a abbreviation column
species_abbrev <- function(full_name, sep_string=". "){
  genus_name <- str_split(full_name, " ")[[1]][1]
  sp_name <- str_split(full_name, " ")[[1]][2]
  genus_abbrev <- str_sub(genus_name, 1, 1)
  combine_name <- paste(genus_abbrev, sp_name, sep = sep_string)
  return (combine_name)
}

symbiosis_tbl[, short_name := mapply(species_abbrev, Species)][]

## read and select data (remove unidentified species as well)
glamap <- read_csv("data/lgm_foram_count_data/GLAMAP.csv") %>% select(!c(`Foram plankt oth [%]`, `Foram plankt [#]`, `Sample comment (Author or institute)`))
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
                                           "G. mentum" = "O. NA"
))


glamap_sps <- unique(glamap$Species)
glamap_sps[!glamap_sps %in% symbiosis_tbl$short_name]

glamap_merged <- merge(glamap, symbiosis_tbl, by.x="Species", by.y = "short_name") %>% select(!c(Species, Species.y))

# aggregate functional group abundance and divided by 100
glamap_merged <- glamap_merged %>% group_by(Campaign, Event, Latitude, Longitude, `Date/Time`, `Depth [m]`, `Elevation [m]`,  Symbiosis, Spinose) %>% 
  summarise_all(.funs = sum, na.rm=T) %>% ungroup() %>% filter(Symbiosis != "Undetermined" & Spinose != 'Undetermined')
glamap_merged <- glamap_merged %>% mutate(Relative_Abundance = Relative_Abundance/100)

# export
write_csv(glamap_merged, "data/lgm_foram_count_data/LGM_GLAMAP_groupped.csv")
