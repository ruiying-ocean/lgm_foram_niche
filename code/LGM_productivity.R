## a R script to pre-process the LGM productivity data

library(data.table)
library(sf)
library(tmap)
library(spData)
library(readxl)
library("tidyverse")

## [] TODO check data converted from image
## [] CANCELLED split comparison of Late Holocene and Pre-Industry

Kohfeld_2013 <- fread("data/LGM_productivity/Kohfeld_etal_2013.csv")
#Radi_2008 <- read_excel("data/LGM_productivity/Radi_etal_2008.xlsx") %>% data.table()
Rae_2020 <- read_excel("data/LGM_productivity/Rae_etal_2020.xlsx") %>% data.table()

#setnames(Radi_2008, "Annual_mean_anomalies(gC/m2)", "Productivity_anomlies")
setnames(Rae_2020, "Sediment_opal_change", "Productivity_anomlies")

# remove NA
Kohfeld_2013 <- Kohfeld_2013[Productivity_anomlies != 99, ]
# from quantitative to qualitative
#Radi_2008 <- Radi_2008[, Productivity_anomlies := ifelse(Productivity_anomlies>0, 1, -1)]

dt_list <- list(Kohfeld_2013, Rae_2020)

for (i in 1:length(dt_list)) {
  select_cols <- c('Core', 'Latitude', 'Longitude', 'Productivity_anomlies')
  dt_list[[i]] <- dt_list[[i]][, ..select_cols]
}
full_data <- rbindlist(dt_list)
string_cols <- c('Latitude', 'Longitude', 'Productivity_anomlies')
full_data[, (string_cols) := lapply(.SD, as.numeric), .SDcols = string_cols]
fwrite(full_data, file="data/LGM_productivity/combined_data.csv")

if (!any(is.na(full_data))) {
  sf_data <- st_as_sf(full_data, coords = c("Longitude", "Latitude"), crs = 4326)
  world <- st_read(system.file("shapes/world.gpkg", package="spData"))
  tm_shape(world)+tm_fill("black")+ tm_shape(sf_data) + tm_dots(col="Productivity_anomlies")
}

#### quantitative NPP
df <- read.csv("data/LGM_productivity/Hernandez-Almeida_etal_2019_raw.csv")
df <- rename(df, Longitude =  Longitude..E., Latitude =  Latitude..N., NPP=NPP.mg.C.m2.day.)
df <- df %>% filter(abs(Latitude) <= 30)
# take average
df <- df %>% group_by(Core, Latitude, Longitude, Time.period) %>% summarise(mean(NPP)) %>% ungroup()
df <- df %>% pivot_wider(id_cols = c("Core", "Latitude", "Longitude"), names_from = "Time.period",values_from = "mean(NPP)")
df <- df %>% mutate(anomaly = LGM - MLH) %>% filter(!is.na(anomaly))
df <- rename(df, anomaly=anomaly(mg_C/m2/day))
write_csv(df, "data/LGM_productivity/Hernandez-Almeida_etal_2019_anomaly.csv")
