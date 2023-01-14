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

## read and calculate relative abundance
climap <- read_csv("data/lgm_foram_count_data/CLIMAP.csv")
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


find_sp <- function(data,symbiosis_tbl_name = 'short_name'){
  sp_list <- data %>% pull(Species) %>% unique()
  unincluded_sp <- sp_list[!sp_list %in% symbiosis_tbl[[symbiosis_tbl_name]]]
  if (length(unincluded_sp) > 0){
    return(unincluded_sp)
  } else{
    print("All species included!")
  }
}
find_sp(climap)

climap_merged <- merge(climap, symbiosis_tbl, by.x="Species", by.y = "short_name") %>% select(!c(Species, Species.y)) %>% select(!`Sample comment`)

# aggregate functional group abundance and divided by 100
climap_merged <- climap_merged %>% group_by(Latitude, Longitude, Event, `Depth [m]`, `Elevation [m]`,  Symbiosis, Spinose) %>% 
  summarise_all(.funs = sum, na.rm=T) %>% ungroup() %>% filter(Symbiosis != "Undetermined" & Spinose != 'Undetermined')

write_csv(climap_merged, "data/lgm_foram_count_data/LGM_CLIMAP_groupped.csv")
