###### Plot site map with chl-a as background
library(tmap)
library(sf)
library(raster)
library(cmocean)
library(tidyverse)

site_info <- read_csv("data/Strack_2022/CoreList_PlanktonicForaminifera.csv")
site_info <- site_info %>% add_row(Site="HU89038-PC8", Lat=33.7, Lon=-57.6, Depth=3500) %>%
  add_row(Site = "KNR140-2-56GGC" , Lat=32.94, Lon=-76.295, Depth=1400) %>%
  add_row(Site = "OMEXII-9K" , Lat=42.33, Lon=-09.465, Depth=1833) %>%
  add_row(Site = "ODP-658" , Lat=20.75, Lon=-34.58, Depth=2263) %>%
  add_row(Site = "MD03-2705" , Lat=18.097, Lon=-21.153, Depth=3085) %>%
  add_row(Site = "GeoB7926-2" , Lat=20.2167, Lon=-18.45, Depth=2500)


extent <- c(-100, -20, 20, 90)

# SeaWifs annual mean chlorophyll a concentration in 2010
chl_a <- brick("data/modern_env/SEASTAR_SEAWIFS_GAC.20100101_20101231.L3m.YR.CHL.chlor_a.9km.nc",
                varname="chlor_a")


p_ocn <-  tm_shape(chl_a, bbox = extent) + tm_raster(breaks = c(0.01, 0.05, 0.1, 1,5, 10, 60),
                                                     style = "order",
                                                     title = "Chl a (mg/m3) 2010 annual mean",
                                                     legend.is.portrait = FALSE, #horizontal legend
                                                     palette = cmocean("haline")(6))
p_ocn <- p_ocn + tm_layout(legend.title.size = .7,legend.outside = F,
                           legend.position = c(0.7, 0.92))

land <- read_sf("data/modern_env/ne_50m_land/ne_50m_land.shp")
p_land <- tm_shape(land, bbox=extent)+ tm_polygons()

site_sf <- st_as_sf(site_info, coords = c("Lon", "Lat"), crs=4326) #WGS84
p_site <- tm_shape(site_sf, bbox=extent) + tm_symbols(shape=18, col="black",size=1)+
  tm_text(text = 'Site', size = .75, auto.placement = T, col="#d1495b") +
  tm_grid(x = seq(-180, 180, by=20), y=seq(-90,90,by=10), 
          lwd = 0.2,col = "gray80")

p <- p_ocn + p_land + p_site

tmap_save(p, "data/Strack_2022/Site_info.png")
