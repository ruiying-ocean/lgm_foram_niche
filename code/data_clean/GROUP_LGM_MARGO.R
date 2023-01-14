library(data.table)
library(dplyr)

## read symbiosis table
source("code/data_clean/read_symbiosis_table.R")

## read and select data (remove unidentified species as well)
north_atlantic <- fread("data/MARGO/MARGO_NorthAtlantic.csv")
south_atlantic <- fread("data/MARGO/MARGO_SouthAtlantic.csv")
pacific <- fread("data/MARGO/MARGO_Pacific_PF.csv")
indopac <- fread("data/MARGO/MARGO_IndoPac.csv")

north_atlantic <- north_atlantic[, c(3:4, 7, 17:53)]
south_atlantic <- south_atlantic[, c(3:4, 7, 17:53)]
pacific <- pacific[, c(3:4, 7, 17:43)]
indopac <- indopac[, c(2:3, 6, 12:48)]

## rename to match symbiosis/spine table
for (DT in list(north_atlantic, south_atlantic, pacific, indopac)){
  setnames(DT, "Globigerinoides sacc", "Globigerinoides sacculifer", skip_absent = T)
  setnames(DT, "Beela digitata", "Globigerinella digitata", skip_absent = T)
  setnames(DT, "Dentagloborotalia anfracta", "Dentigloborotalia anfracta", skip_absent = T)
  setnames(DT, "Neogloboquadrina pachyderma L", "Neogloboquadrina pachyderma", skip_absent = T)
  setnames(DT, "Neogloboquadrina pachyderma R", "Neogloboquadrina incompta", skip_absent = T)
  setnames(DT, "Globorotalia menardii flexuosa", "Globorotalia flexuosa", skip_absent = T)
  setnames(DT, "Globorotalia crassula","Globorotalia crassa", skip_absent = T)
}
remove(DT)

north_atlantic <- north_atlantic[!is.na(Longitude) & !is.na(Latitude), ]
south_atlantic <- south_atlantic[!is.na(Longitude) & !is.na(Latitude), ]
pacific <- pacific[!is.na(Longitude) & !is.na(Latitude), ]
indopac <- indopac[!is.na(Longitude) & !is.na(Latitude), ]

## sum across species columns and check is count or percentage data
north_atlantic[, `:=`(total_foram = rowSums(.SD, na.rm=T)), .SDcols=4:length(north_atlantic)]
south_atlantic[, `:=`(total_foram = rowSums(.SD, na.rm=T)), .SDcols=4:length(south_atlantic)]
pacific[, `:=`(total_foram = rowSums(.SD, na.rm=T)), .SDcols=4:length(pacific)]
indopac[, `:=`(total_foram = rowSums(.SD, na.rm=T)), .SDcols=4:length(indopac)]

# if count data, then all columns %% 1 == 0
indopac_israw  <- indopac %>%
  mutate(across(4:41, ~.%%1)) %>%
  mutate(is_raw = select(., 4:41) %>%
           rowSums(na.rm = TRUE)) %>%
  mutate(is_raw = is_raw == 0) %>%
  select(is_raw)

indopac[, is_raw := indopac_israw$is_raw]
indopac_raw <- indopac[is_raw == TRUE,]
indopac_perc <- indopac[is_raw != TRUE, ]
remove(indopac_israw)

## convert raw data to relative abundance
indopac_raw <- indopac_raw |> mutate(across(4:40, ~./total_foram))

indopac_raw <- indopac_raw[, !c("total_foram", "is_raw", "Other ID")]
indopac_perc <- indopac_perc[, !c("total_foram", "is_raw", "Other ID")]

indopac <- indopac[, !c("total_foram", "is_raw", "Other ID")]
pacific <- pacific[, !("total_foram")]
south_atlantic <- south_atlantic[, !c("total_foram", "Other ID", "Other UnID")]
north_atlantic <- north_atlantic[, !c("total_foram", "Other ID", "Other UnID")]

## check all species are included
which(!(names(indopac)[4:39] %in% symbiosis_tbl$Species))
which(!(names(north_atlantic)[4:38] %in% symbiosis_tbl$Species))
which(!(names(pacific)[4:30] %in% symbiosis_tbl$Species))
which(!(names(south_atlantic)[4:38] %in% symbiosis_tbl$Species))

## covert them to long format
pacific <- melt.data.table(pacific, id.vars = c("Latitude", "Longitude", "Sample depth - upper (m)"), variable.name = "Species")
north_atlantic <- melt.data.table(north_atlantic, id.vars = c("Latitude", "Longitude", "Sample depth - upper (m)"), variable.name = "Species")
south_atlantic <- melt.data.table(south_atlantic, id.vars = c("Latitude", "Longitude", "Sample depth - upper (m)"),  variable.name = "Species")
indopac_perc <- melt.data.table(indopac_perc, id.vars = c("Latitude", "Longitude", "Sample depth - upper (m)"),  variable.name = "Species")
indopac_raw <- melt.data.table(indopac_raw, id.vars = c("Latitude", "Longitude", "Sample depth - upper (m)" ), variable.name = "Species", )

## join symbiosis data
pacific_merged <- merge(pacific, symbiosis_tbl, by="Species", all.x = TRUE)
south_atlantic_merged <- merge(south_atlantic, symbiosis_tbl, by="Species", all.x = TRUE)
north_atlantic_merged <- merge(north_atlantic, symbiosis_tbl, by="Species", all.x = TRUE)
indopac_raw_merged <- merge(indopac_raw, symbiosis_tbl, by="Species", all.x = TRUE)
indopac_perc_merged <- merge(indopac_perc, symbiosis_tbl, by="Species", all.x = TRUE)

pacific_merged <- pacific_merged[!is.na(value),][,!"Species"]
south_atlantic_merged <- south_atlantic_merged[!is.na(value),][,!"Species"]
north_atlantic_merged <- north_atlantic_merged[!is.na(value),][,!"Species"]
indopac_raw_merged <- indopac_raw_merged[!is.na(value),][,!"Species"]
indopac_perc_merged <- indopac_perc_merged[!is.na(value),][,!"Species"]

## aggregate (note indopac_raw doesn't need to be divided by 100)
pacific_merged <- pacific_merged[, .(Prop = sum(value)/100), by = .(Latitude, Longitude,`Sample depth - upper (m)`, Symbiosis,Spinose)]
indopac_raw_merged <- indopac_raw_merged[, .(Prop = sum(value)), by = .(Latitude, Longitude,`Sample depth - upper (m)`, Symbiosis,Spinose)]
indopac_perc_merged <- indopac_perc_merged[, .(Prop = sum(value)/100), by = .(Latitude, Longitude,`Sample depth - upper (m)`, Symbiosis,Spinose)]
south_atlantic_merged <- south_atlantic_merged[, .(Prop = sum(value)/100), by = .(Latitude, Longitude,`Sample depth - upper (m)`, Symbiosis,Spinose)]
north_atlantic_merged <- north_atlantic_merged[, .(Prop = sum(value)/100), by = .(Latitude, Longitude,`Sample depth - upper (m)`, Symbiosis,Spinose)]

## merge oceans
margo_data <- rbindlist(list(north_atlantic_merged, south_atlantic_merged,
                             pacific_merged, indopac_perc_merged, indopac_raw_merged))

margo_data <- margo_data[Longitude <= 180, ]
fwrite(margo_data, "data/MARGO/MARGO_grouped.csv")
