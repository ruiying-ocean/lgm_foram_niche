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
epilog <- read_csv("data/lgm_foram_count_data/LGM_EPILOG.csv") %>% select(-`Foram plankt [#]`)
epilog <- epilog %>% pivot_longer(cols=c(`G. bulloides [%]`:`T. iota [%]`), names_to = "Species", values_to="Relative_Abundance")
epilog <- epilog %>%  mutate_at("Species", str_replace, " \\[%\\]", "")

## Rename species to match symbiosis table
### Get the different names
#epilog_sps <- unique(epilog$Species)
#epilog_sps[!epilog_sps %in% symbiosis_tbl$short_name]

epilog <- epilog %>% mutate(Species=recode(Species, 
                                 "N. pachyderma d" =  "N. pachyderma",
                                 "N. pachyderma s" = "N. incompta",
                                 "G. ruber w" = "G. ruber",
                                 "G. ruber p" = "G. ruber",
                                 "G. ruber hsp" = "G. ruber",
                                 "G. menardii flexuosa" = "G. menardii",
                                 "G. rubescens white" = "G. rubescens",
                                 "G. rubescens pink" = "G. rubescens",
                                 "G. truncatulinoides s" = "G. truncatulinoides",
                                 "G. truncatulinoides d" = "G. truncatulinoides",
                                 "G. sacculifer wo sac" = "G. sacculifer",
                                 "G. sacculifer sac" = "G. sacculifer",
                                 "G. quinqueloba"= "T. quinqueloba",
                                 "G. bradyi" = "G. uvula",
                                 "P/D int"="N. dutertrei",
                                 ))
epilog_merged <- merge(epilog, symbiosis_tbl, by.x="Species", by.y = "short_name") %>% select(!c(Species, Species.y))

# aggregate functional group abundance and divided by 100
epilog_merged <- epilog_merged %>% group_by(Latitude, Longitude, Event, `Depth [m]`, `Elevation [m]`,  Symbiosis, Spinose) %>% 
  summarise_all(.funs = sum, na.rm=T) %>% ungroup() %>% filter(Symbiosis != "Undetermined" & Spinose != 'Undetermined')
epilog_merged <- epilog_merged %>% mutate(Relative_Abundance = Relative_Abundance/100)

# export
write_csv(epilog_merged, "data/lgm_foram_count_data/LGM_EPILOG_groupped.csv")
