# No-analogue assemblages to LGM (Strack et al., 2022, NEE)
# Anne Strack, last edit 28 July 2022

# Data citation:
# Kucera, M., Rosell-Melé, A., Schneider, R., Waelbroeck, C. & Weinelt, M. Multiproxy approach for the reconstruction 
# of the glacial ocean surface (MARGO). Quat. Sci. Rev. 24, 813-819, doi:10.1016/j.quascirev.2004.07.017 (2005).
# Kucera, M. et al. Reconstruction of sea-surface temperatures from assemblages of planktonic foraminifera: multi-
# technique approach based on geographically constrained calibration data sets and its application to glacial Atlantic 
# and Pacific Oceans. Quat. Sci. Rev. 24, 951-998, doi:10.1016/j.quascirev.2004.07.014 (2005).

# Data link: 
# Mediterranean: https://doi.pangaea.de/10.1594/PANGAEA.227306
# Atlantic:      https://doi.pangaea.de/10.1594/PANGAEA.227329
# Pacific:       https://doi.pangaea.de/10.1594/PANGAEA.227327

# This script calculates the compositional dissimilarity to the nearest LGM sample (results shown in Figure 4) to 
# evaluate the existence of no-analogue assemblages. It also produces the raw plots of Figure 4 as well as 
# Extended Data Fig. 3 and 4. Final figure preparation was conducted in CorelDraw.


remove(list = ls(all.names = TRUE)) # clear environment including hidden objects

# --- load packages
library(tidyverse)
library(janitor) # used for cleaning tables
library(vegan) # Community Ecology Package: provides tools for descriptive community ecology
library(codyn) # Community Dynamics Metrics Package
library(viridis)
library(gridExtra) # used for grid.arrange
library(cowplot) # 1.0.0  
library(matrixStats) # used for i.e. rowMedians
library(rioja)
library(data.table) # used for setcolorder
library(raster) # used for gridding data (if this package is used select is masked --> use dplyr::select to avoid error message)
library(directlabels)

# ------------------------------------------------------------------------------------------------------------------------
# --- A. set paths for source files and variables
# ------------------------------------------------------------------------------------------------------------------------
core.list.data = "CoreList_PlanktonicForaminifera.csv"
reference.list = "ReferenceList_PlanktonicForaminifera.csv"
full.data = "FullDataTable_PF_harmonized.txt"
plankton.group = "planktonic foraminifera"

lgm.data = "MARGO_LGM_dataset.csv"

# --- set variables 
min_age <- 0 # minimum age used for analyses
max_age <- 24 # maximum age used for analyses
seq.plots <- seq(0, 24, 3) # set x-axis ticks for plots

# --- set string with possible NA spellings
na_strings <- c("NA", "N A", "N / A", "N/A", "N/ A", "#NA", "#N/A", "nan", "NaN") # string with possible NA spellings (used for data cleaning)


# ------------------------------------------------------------------------------------------------------------------------
# --- B. load and clean core list
# ------------------------------------------------------------------------------------------------------------------------
df_coreList <- read.csv(core.list.data, na = na_strings) %>% 
  clean_names() %>% # uses janitor to clean all names (all lower case with underscore)
  remove_empty(which = c("rows", "cols")) %>%  # removes all empty rows and columns (sometimes produced by Excel)
  filter(plankton == plankton.group) %>% # selects only sites that are used in analysis
  mutate("ID" = str_c(core, plankton, sep = "_")) %>% # add ID column do core list (combination of site & plankton group)
  dplyr::select(ID, everything()) # set ID column as first column


# -----------------------------------------------------------------------------------------------------------------------
# --- C. Functions used in this script
# -----------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------
# --- D. read reference list of PF names
# ------------------------------------------------------------------------------------------------------------------------
df_reference <- read.csv(reference.list) %>%
  remove_empty(which = c("rows", "cols"))  # removes all empty rows and columns (sometimes produced by Excel)
reference <- setNames(rep(NA, length(df_reference$Species)), df_reference$Species) # produces numeric vector of reference species


# ------------------------------------------------------------------------------------------------------------------------
# --- E. load full data table and adjust col names
# ------------------------------------------------------------------------------------------------------------------------
# --- read fully denormalized data table
df_full <- read_tsv(full.data) %>%
  filter(Age_kaBP >= min_age & Age_kaBP <= max_age & Plankton == "planktonic foraminifera") # only PF used, only last 24 ka

# --- produce wider table
df_filtered_wider <- df_full %>%
  pivot_wider(names_from = "Species", values_from = "Rel_abundance") # makes table wider again

# --- Which sites have combined G. ruber, combined G.menardii/tumida?
comb.ruber <- unique(df_full$ID[df_full$Species == "Globigerinoides_ruber_ruber_and_Globigerinoides_ruber_albus"])
comb.men.tum <- unique(df_full$ID[df_full$Species == "Globorotalia_menardii_and_Globorotalia_tumida"])
# --- deleting combined Globorotalia_menardii_and_Globorotalia_tumida column (only site that has this column is MD99-2284 and it only contains zeros)
df_filtered_wider <- df_filtered_wider %>%
  dplyr::select(-Globorotalia_menardii_and_Globorotalia_tumida)
# --- summing up G.ruber ruber and albus and remove individual columns (as some sites only contain combined G. ruber)
df_filtered_wider <- df_filtered_wider %>%
  mutate(Globigerinoides_ruber_ruber_and_Globigerinoides_ruber_albus = 
           ifelse(is.na(df_filtered_wider$Globigerinoides_ruber_ruber_and_Globigerinoides_ruber_albus), 
                  Globigerinoides_ruber_ruber + Globigerinoides_ruber_albus, 
                  Globigerinoides_ruber_ruber_and_Globigerinoides_ruber_albus)) %>%
  dplyr::select(-Globigerinoides_ruber_ruber, -Globigerinoides_ruber_albus)

# --- adjusting col names of my data set to reference list
# adding missing columns
df_filtered_wider <- df_filtered_wider %>%
  add_column(!!!reference[!names(reference) %in% names(.)]) # adds non-existing cols from ref to df with NA
# reordering cols according to reference list
df_filtered_wider_reordered <- setcolorder(df_filtered_wider[-(1:5)], names(reference)) %>% # reorders col names using reference
  mutate_all(funs(./100)) %>%
  add_column(ID = df_filtered_wider$ID, Site = df_filtered_wider$Site, Plankton = df_filtered_wider$Plankton, Depth_m = df_filtered_wider$Depth_m, Age_kaBP = df_filtered_wider$Age_kaBP) %>% # adds core and lat data
  dplyr::select(ID, Site, Plankton, Depth_m, Age_kaBP, everything())  

df_filtered_wider_reordered[is.na(df_filtered_wider_reordered)] <- 0 # replaces NA with 0


# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
#### --- START: No-analogue community analysis
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------
# --- load and adjust MARGO LGM data set
# ------------------------------------------------------------------------------------------------------------------------
# download MARGO data set from https://doi.pangaea.de/10.1594/PANGAEA.227306 (Mediterranean),
#                              https://doi.pangaea.de/10.1594/PANGAEA.227329 (Atlantic), and
#                              https://doi.pangaea.de/10.1594/PANGAEA.227327 (Pacific)
# save all subsets manually in one .csv file named "MARGO_LGM_dataset.csv" or adjust file name

df_lgm <- read.csv(lgm.data, na = na_strings, stringsAsFactors = FALSE) %>%
  clean_names() %>% # uses janitor to clean all names (all lower case with underscore)
  remove_empty(which = c("rows", "cols")) %>%  # removes all empty rows and columns (sometimes produced by Excel)
  filter(ocean == "North Atlantic" | ocean == "South Atlantic" | ocean == "Mediterranean", lat >= -6)  %>% # only uses Atlantic data with same latitudinal extent than my data
  dplyr::select(-c(count, type, sum)) %>%
  rename(Globigerinoides_ruber_ruber_and_Globigerinoides_ruber_albus = "globigerinoides_ruber_ruber_and_globigerinoides_ruber_albus") #%>% # rename this column as next line capitalizes only first letter of col name
names(df_lgm) <- paste(toupper(substring(names(df_lgm), 1, 1)), substring(names(df_lgm), 2), sep = "") # capitalizes first letter of column names

# --- summing up G.ruber ruber and albus and remove individual columns
df_lgm <- df_lgm %>%
  mutate(Globigerinoides_ruber_ruber_and_Globigerinoides_ruber_albus =
           ifelse(is.na(df_lgm$Globigerinoides_ruber_ruber_and_Globigerinoides_ruber_albus),
                  Globigerinoides_ruber_ruber + Globigerinoides_ruber_albus,
                  Globigerinoides_ruber_ruber_and_Globigerinoides_ruber_albus)) %>%
  dplyr::select(-Globigerinoides_ruber_ruber, -Globigerinoides_ruber_albus)

# --- adjusting col names of MARGO data set to reference list
# adding missing columns
df_lgm <- df_lgm %>%
  dplyr::select(-Unidentified) %>% # removed unidentified species, as I also removed them from my data set --> recalculation to 100 later on.
  add_column(!!!reference[!names(reference) %in% names(.)]) # adds non-existing cols from ref to df with NA

# reordering cols according to reference list
df_lgm_reordered <- setcolorder(df_lgm[-(1:15)], names(reference)) %>% # reorders col names using reference
  add_column(Core = df_lgm$Core, Lat = df_lgm$Lat, Long = df_lgm$Long, Water_depth_m = df_lgm$Water_depth_m,
             Depth_upper_m = df_lgm$Sample_depth_upper_m, Depth_lower_m = df_lgm$Sample_depth_lower_m,
             Age_kaBP = df_lgm$Calendar_age_estimate_cal_ky_bp) %>% # adds core and lat data
  dplyr::select(Core, Lat, Long, Water_depth_m, Depth_upper_m, Depth_lower_m, Age_kaBP, everything()) %>%
  transform(Age_kaBP = as.numeric(Age_kaBP))
df_lgm_reordered[-(1:7)][is.na(df_lgm_reordered[-(1:7)])] <- 0 # replaces NA with 0 (but not metadata columns)

# --- percent calculation and recalculation to sum up to 1 (removed unidentified col beforehand)
df_lgm_reordered <- modify_at(df_lgm_reordered, names(df_lgm_reordered[-(1:7)]), ~./rowSums(df_lgm_reordered[-(1:7)]))


# ------------------------------------------------------------------------------------------------------------------------
# --- update MARGO LGM data set
# ------------------------------------------------------------------------------------------------------------------------
# sites from this compilation that are not yet included in the LGM data set
CoresNotInMARGO <- setdiff(df_filtered_wider_reordered$Site, df_lgm_reordered$Core)
CoresInMARGO <- setdiff(df_filtered_wider_reordered$Site, CoresNotInMARGO)

### --- produces df that only contains LGM samples of this study that are not included in MARGO (n = 194)
df_CoresNotInMargo <- df_filtered_wider_reordered %>% 
  filter(Site %in% CoresNotInMARGO) %>% # only show cores of this study that are not included in MARGO
  filter(Age_kaBP >= 19 & Age_kaBP <= 23) %>% # select LGM samples (interval defined by EPILOG Workshop (Mix et al., 2001))
  left_join(df_coreList[, c("ID", "lat", "lon", "depth")], by = "ID") %>% # looks up and adds coordinates to df
  dplyr::select(-c(ID, Plankton)) %>%
  mutate(Depth_lower_m = NA) %>%
  dplyr::select(Site, lat, lon, depth, Depth_m, Depth_lower_m, Age_kaBP, everything()) %>%
  rename("Core" = Site, "Lat" = lat, "Long" = lon, "Depth_upper_m" = Depth_m, "Water_depth_m" = depth)
### --- produce df that contains ALL LGM samples (19-23 ka) of this study no matter if part of MARGO or not (n = 288)
df_myCoresLGM <- df_filtered_wider_reordered %>% # only show cores that are included in LGM
  filter(Site %in% unique(df_filtered_wider_reordered$Site)) %>%
  filter(Age_kaBP >= 19 & Age_kaBP <= 23) %>% # select LGM samples (interval defined by EPILOG Workshop (Mix et al., 2001))
  left_join(df_coreList[, c("ID", "lat", "lon", "depth")], by = "ID") %>% # looks up and adds coordinates to df
  dplyr::select(-c(ID, Plankton)) %>%
  mutate(Depth_lower_m = NA) %>%
  dplyr::select(Site, lat, lon, depth, Depth_m, Depth_lower_m, Age_kaBP, everything()) %>%
  rename("Core" = Site, "Lat" = lat, "Long" = lon, "Depth_upper_m" = Depth_m, "Water_depth_m" = depth)
### --- removes LGM samples from sites used in this study that are already in MARGO and adds/updates my LGM samples
df_lgm_reordered <- df_lgm_reordered[!(df_lgm_reordered$Core %in% CoresInMARGO),] # remove LGM samples from sites that are in this data set
df_lgm_reordered <- rbind(df_lgm_reordered, df_myCoresLGM) # add and/or update LGM samples from my data


# ------------------------------------------------------------------------------------------------------------------------
# --- produce overview map of LGM data set (Extended Data Fig. 1)
# ------------------------------------------------------------------------------------------------------------------------
map <- map_data("world")

map_lgm <- ggplot() +
  geom_polygon(data = map, aes(long, lat, group = group), fill = "grey95", colour = "grey95") +
  coord_fixed(xlim = c(-80, 40), ylim = c(-5, 83)) +
  labs(x = "Longitude", y = "Latitude", title = "LGM dataset (MARGO)") +
  theme(panel.background = element_rect(fill = "grey80"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "right",
        legend.direction ="vertical",
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(colour = "black", size = 10),
        axis.text.y = element_text(colour = "black", size = 10),
        plot.title = element_blank(),
        legend.text = element_text(size = 10)) +
  geom_point(data = unique(dplyr::select(df_lgm_reordered, "Core", "Lat", "Long")), aes(x = Long, y = Lat), size = 1) +
  geom_point(data = unique(dplyr::select(df_CoresNotInMargo, "Core", "Lat", "Long")), aes(x = Long, y = Lat), size = 1, colour = "red3") +
  scale_x_continuous(n.breaks = 6) # expand is used to change the padding around the data (in percent for continious data)
map_lgm

ggsave(filename = paste0("ED_Fig3_OverviewMap_LGMdataset.png"), plot = map_lgm, units = "cm", height = 20, width = 30, dpi = 300)


# ------------------------------------------------------------------------------------------------------------------------
# --- calculate minimal distance of each sample in this study to nearest LGM sample
# ------------------------------------------------------------------------------------------------------------------------
# - calculate M-H index with vegan
df_filtered_lgm <- bind_rows(df_filtered_wider_reordered[-(1:5)], df_lgm_reordered[-(1:7)]) # combine my data with lgm data to one df

horn <- vegdist(x = df_filtered_lgm, method = "horn") # calculates M-H index
df_horn <- as_tibble(as(horn, "matrix")) # make tibble
dist_horn <- apply(df_horn[1:nrow(df_filtered_wider_reordered), -(1:nrow(df_filtered_wider_reordered))], 1, min) # subset df_Horn and look for minimal distance

# --- join distance and lat/long information to data
df_filtered_wider_reordered <- df_filtered_wider_reordered %>%
  mutate(Dist_horn = dist_horn) %>% # adds min distances to df
  left_join(df_coreList[, c("ID", "lat", "lon")], by = "ID") # looks up and adds coordinates to df


# ------------------------------------------------------------------------------------------------------------------------
#### --- plot results for nearest analogue to LGM data set (Figure 4)
# ------------------------------------------------------------------------------------------------------------------------
df_raster <- df_filtered_wider_reordered %>%
  dplyr::select(x = Age_kaBP, y = lat, z = Dist_horn)

e <- extent(c(0, 24, -10, 75))
r <- raster(e, ncol=24, nrow=34)

r <- rasterize(df_raster[,1:2], r, df_raster[,3], fun = mean)

raster <- as.data.frame(r, xy = TRUE)


plot.hovm <- ggplot(data = df_filtered_wider_reordered) +
  labs(x = "Age [ka]", y = "Latitude"#,
       #title = "Hovmoller diagram: Morisita-Horn to nearest LGM sample"
       ) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "right",
        legend.direction ="vertical",
        legend.background = element_blank(),
        legend.title = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10),
        axis.title.y = element_text(colour = "black", size = 12),
        axis.title.x = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 10),
        axis.text.y = element_text(colour = "black", size = 10),
        plot.title = element_text(size = 12)) +
  scale_x_continuous(breaks = seq.plots, limits = c(min_age, max_age), expand = c(0.005,0.005)) + # expand is used to change the padding around the data (in percent for continious data)
  scale_y_continuous(limits = c(-10, 75), expand = c(0.001,0.001)) + # expand is used to change the padding around the data (in percent for continious data)
  geom_raster(data = raster, aes(x = x, y = y, fill = layer),
              interpolate = FALSE) +
  geom_point(aes(x = Age_kaBP, y = lat), colour = "grey80", size = 1) +
  scale_fill_gradient2(low = "gray98", mid = "gray90", midpoint = 0.06, high = "midnightblue", na.value = "white", name = "Compositional\ndissimilarity") # use this option for raw data maps
plot.hovm

ggsave(filename = "Figure4_NoAnalogues.png", plot = plot.hovm, units = "cm", height = 10, width = 18, dpi = 300)

# ------------------------------------------------------------------------------------------------------------------------
# --- calculate distribution of pair-wise dissimilarities within LGM database (Extended Data Figure 2)
# ------------------------------------------------------------------------------------------------------------------------

# --- calculate M-H index for LGM dataset
horn_lgm <- vegdist(x = df_lgm_reordered[-(1:7)], method = "horn", upper = FALSE) # calculates M-H
df_horn_lgm <- as_tibble(as(horn_lgm, "matrix")) # make tibble

dist_horn_lgm <- df_horn_lgm %>%
  mutate(nearest1 = apply(df_horn_lgm, 1, function(x) sort(x)[2]), # distance to nearest non-self LGM sample
         nearest2 = apply(df_horn_lgm, 1, function(x) sort(x)[3]), # distance to 2nd nearest non-self LGM sample
         nearest3 = apply(df_horn_lgm, 1, function(x) sort(x)[4])) # distance to 3rd nearest non-self LGM sample

# --- calculate 95 and 99 percentiles
quan1 <- quantile(dist_horn_lgm$nearest1, c(0.95,0.99))
quan2 <- quantile(dist_horn_lgm$nearest2, c(0.95,0.99))
quan3 <- quantile(dist_horn_lgm$nearest3, c(0.95,0.99))
quan <- quantile(df_filtered_wider_reordered$Dist_horn, c(0.95,0.99))

# ------------------------------------------------------------------------------------------------------------------------
# --- plot results for distribution of pair-wise dissimilarities within LGM database (Extended Data Figure 2)
# ------------------------------------------------------------------------------------------------------------------------

# --- histograms within LGM dataset: distance to nearest non-self LGM sample (Extended Data Figure 2)
plot.hist1 <- ggplot(data = dist_horn_lgm) +
  geom_histogram(aes(x = nearest1),
                 binwidth = 0.01, boundary = 0,
                 fill = "grey70", colour = NA) +
  geom_vline(aes(xintercept = quan1[1]),
             color = "black", linetype = 3, size = 1) +
  annotate("text", x = quan1[1], y = 500, label = paste0("95 percentile: ", sprintf("%.3f", round(quan1[1], digits = 3)),"\n"),
           colour="black", angle = 90) +
  geom_vline(aes(xintercept = quan1[2]),
             color = "black", linetype = 2, size = 1) +
  annotate("text", x = quan1[2], y = 500, label = paste0("99 percentile: ", sprintf("%.3f", round(quan1[2], digits = 3)),"\n"),
           colour="black", angle = 90) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "right",
        legend.direction ="vertical",
        legend.background = element_blank(),
        legend.title = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10),
        axis.title.y = element_text(colour = "black", size = 12),
        axis.title.x = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 10),
        axis.text.y = element_text(colour = "black", size = 10),
        plot.title = element_text(size = 12)) +
  scale_x_continuous(breaks = seq(0,1, 0.1), limits = c(0,0.71), expand = c(0.01,0.01)) + # expand is used to change the padding around the data (in percent for continious data)
  labs(x = "Dissimilarity", y = "Count",
       title = "Distances to nearest non-self sample within LGM database")
plot.hist1

plot.hist2 <- ggplot(data = dist_horn_lgm) +
  geom_histogram(aes(x = nearest2),
                 binwidth = 0.01, boundary = 0,
                 fill = "grey70", colour = NA) +
  geom_vline(aes(xintercept = quan2[1]),
             color = "black", linetype = 3, size = 1) +
  annotate("text", x = quan2[1], y = 400, label = paste0("95 percentile: ", sprintf("%.3f", round(quan2[1], digits = 3)),"\n"),
           colour="black", angle = 90) +
  geom_vline(aes(xintercept = quan2[2]),
             color = "black", linetype = 2, size = 1) +
  annotate("text", x = quan2[2], y = 400, label = paste0("99 percentile: ", sprintf("%.3f", round(quan2[2], digits = 3)),"\n"),
           colour="black", angle = 90) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "right",
        legend.direction ="vertical",
        legend.background = element_blank(),
        legend.title = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10),
        axis.title.y = element_text(colour = "black", size = 12),
        axis.title.x = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 10),
        axis.text.y = element_text(colour = "black", size = 10),
        plot.title = element_text(size = 12)) +
  scale_x_continuous(breaks = seq(0,1, 0.1), limits = c(0,0.71), expand = c(0.01,0.01)) + # expand is used to change the padding around the data (in percent for continious data)
  labs(x = "Dissimilarity", y = "Count",
       title = "Distances to 2nd nearest non-self sample within LGM database")
plot.hist2

plot.hist3 <- ggplot(data = dist_horn_lgm) +
  geom_histogram(aes(x = nearest3),
                 binwidth = 0.01, boundary = 0,
                 fill = "grey70", colour = NA) +
  geom_vline(aes(xintercept = quan3[1]),
             color = "black", linetype = 3, size = 1) +
  annotate("text", x = quan3[1], y = 350, label = paste0("95 percentile: ", sprintf("%.3f", round(quan3[1], digits = 3)),"\n"),
           colour="black", angle = 90) +
  geom_vline(aes(xintercept = quan3[2]),
             color = "black", linetype = 2, size = 1) +
  annotate("text", x = quan3[2], y = 350, label = paste0("99 percentile: ", sprintf("%.3f", round(quan3[2], digits = 3)),"\n"),
           colour="black", angle = 90) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "right",
        legend.direction ="vertical",
        legend.background = element_blank(),
        legend.title = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 10),
        axis.title.y = element_text(colour = "black", size = 12),
        axis.title.x = element_text(colour = "black", size = 12),
        axis.text.x = element_text(colour = "black", size = 10),
        axis.text.y = element_text(colour = "black", size = 10),
        plot.title = element_text(size = 12)) +
  scale_x_continuous(breaks = seq(0,1, 0.1), limits = c(0,0.71), expand = c(0.01,0.01)) + # expand is used to change the padding around the data (in percent for continious data)
  labs(x = "Dissimilarity", y = "Count",
       title = "Distances to 3rd nearest non-self sample within LGM database")
plot.hist3

p <- cowplot::plot_grid(plot.hist1, plot.hist2, plot.hist3, align = "vh", ncol = 1, axis = "tbl")

ggsave(filename = "ED_Fig4_Histogram_Threshold_NoAnalogues.png", plot = p, units = "cm", height = 25, width = 20, dpi = 300)

# # --- histogram distance to nearest LGM analogue (This is the distribution histogram of Figure 4, not shown in manuscript.)
# plot.hist <- ggplot(data = df_filtered_wider_reordered) +
#   geom_histogram(aes(x = Dist_horn),
#                  binwidth = 0.01, boundary = 0,
#                  fill = "grey70", colour = NA) +
#   geom_vline(aes(xintercept = quan[1]),
#              color = "black", linetype = 3, size = 1) +
#   annotate("text", x = quan[1], y = 500, label = paste0("95 percentile: ", sprintf("%.3f", round(quan[1], digits = 3)),"\n"),
#            colour="black", angle = 90) +
#   geom_vline(aes(xintercept = quan[2]),
#              color = "black", linetype = 2, size = 1) +
#   annotate("text", x = quan[2], y = 500, label = paste0("99 percentile: ", sprintf("%.3f", round(quan[2], digits = 3)),"\n"),
#            colour="black", angle = 90) +
#   theme(panel.background = element_rect(fill = "white"),
#         panel.grid = element_blank(),
#         panel.border = element_rect(fill = NA, colour = "black"),
#         legend.position = "right",
#         legend.direction ="vertical",
#         legend.background = element_blank(),
#         legend.title = element_text(colour = "black", size = 10),
#         legend.text = element_text(size = 10),
#         axis.title.y = element_text(colour = "black", size = 12),
#         axis.title.x = element_text(colour = "black", size = 12),
#         axis.text.x = element_text(colour = "black", size = 10),
#         axis.text.y = element_text(colour = "black", size = 10),
#         plot.title = element_text(size = 12)) +
#   scale_x_continuous(breaks = seq(0,1, 0.05), limits = c(0,0.3), expand = c(0.01,0.01)) + # expand is used to change the padding around the data (in percent for continious data)
#   labs(x = "Dissimilarity", y = "Count",
#        title = "Distances to nearest LGM analogue")
# plot.hist
