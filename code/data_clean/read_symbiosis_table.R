library(data.table)

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


find_sp <- function(data,symbiosis_tbl_name = 'short_name'){
  sp_list <- data %>% pull(Species) %>% unique()
  unincluded_sp <- sp_list[!sp_list %in% symbiosis_tbl[[symbiosis_tbl_name]]]
  if (length(unincluded_sp) > 0){
    return(unincluded_sp)
  } else{
    print("All species included!")
  }
}
