# Principal component analysis (PCA) of PF assemblage data (Strack et al., 2022, NEE)
# Anne Strack, last edit 28 July 2022

# Data citation:
# 1) Osman, M. B. et al. Globally resolved surface temperatures since the Last Glacial Maximum. 
#    Nature 599, 239-244, doi:10.1038/s41586-021-03984-4 (2021).
# 2) Locarnini, R. A. et al. World Ocean Atlas 2018, Volume 1: Temperature. A. Mishonov, Technical Editor. 
#    NOAA Atlas NESDIS 81, 52 (2019).

# Data link: 
# 1) https://www.ncei.noaa.gov/pub/data/paleo/reconstructions/osman2021/
# 2) https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/bin/woa18.pl


# This script runs the PCA on the assemblage data of the individual sites (results shown in Figure 2)
# and also on the whole dissimilarity matrix (results shown in Figure 1b). This script also produces the
# raw plots of Figures 1 and 2 of the manuscript. Final figure preparation was conducted in CorelDraw.


remove(list = ls(all.names = TRUE)) # clear environment including hidden objects

# --- load packages
library(tidyverse)
library(janitor) # used for cleaning tables
library(FactoMineR) # used for PCA analysis
library(factoextra) # used for PCA data visualization
library(vegan) # Community Ecology Package: provides tools for descriptive community ecology
library(viridis)
library(broom) # used for "augment" function (see loess model)
library(ncdf4) # package for netcdf manipulation
library(raster) # used for gridding data (if this package is used select is masked --> use dplyr::select to avoid error message)
library(data.table) # used for setcolorder

# ------------------------------------------------------------------------------------------------------------------------
# --- A. set paths for source files and variables
# ------------------------------------------------------------------------------------------------------------------------
core.list.data = "CoreList_PlanktonicForaminifera.csv"
reference.list = "ReferenceList_PlanktonicForaminifera.csv"
full.data = "FullDataTable_PF_harmonized.txt"
plankton.group = "planktonic foraminifera"
temp.data = "GlobalMeanSST_Osman_etal_2021.csv"
WOA18.data = "WOA18_SST_decav_t00mn04.csv"

# --- set variables 
min_age <- 0 # minimum age used for analyses
max_age <- 24 # maximum age used for analyses
seq.plots <- seq(0, 24, 3) # set x-axis ticks for plots
seq.interpol <- seq(0,24,0.5) # time steps used for interpolating all PC1, so they have same time steps (500 year bin)
t.PC1 <- tibble(Age_kaBP = seq.interpol) # produces tibble with time steps set by "sequence"

# --- set string with possible NA spellings
na_strings <- c("NA", "N A", "N / A", "N/A", "N/ A", "#NA", "#N/A", "nan", "NaN") # string with possible NA spellings (used for data cleaning)

# ------------------------------------------------------------------------------------------------------------------------
# --- B. load and clean core list
# ------------------------------------------------------------------------------------------------------------------------
df_coreList <- read.csv(core.list.data, na = na_strings) %>% 
  clean_names() %>% # uses janitor to clean all names (all lower case with underscore)
  remove_empty(which = c("rows", "cols")) %>%  # removes all empty rows and columns (sometimes produced by Excel)
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
# --- START WITH ANALYSIS
# ------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
# --- PCA analysis (shown in Figure 2)
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
id <- unique(df_filtered_wider_reordered$ID) # lists all IDs in a vetor
site <- unique(df_filtered_wider_reordered$Site) # lists all sites in a vetor

# --- write function to run PCA on site-specific assemblage data
PCA.analysis <- function(id) {
  #id <- "161-977A_planktonic foraminifera" # [uncomment if you want to test single site]
  # --- filter data by ID
  df_filtered_wider_reordered <- filter(df_filtered_wider_reordered, ID == id)

  # --- PCA analysis (PCA() automatically standardizes data, no need to do it manually)
  res.pca <- PCA(df_filtered_wider_reordered[-(1:5)], ncp = 3, graph = FALSE) # PCA on dataset is calculated

  eig.val <- get_eigenvalue(res.pca) # eigenvalues of each PC (eig.val[1, 2] gives % of var of PC1)

  #scree.plot <- fviz_eig(res.pca, addlabels = TRUE) # shows % of variance explained for each PC
  #corr.plot <- corrplot(res.pca$var$contrib, is.corr = FALSE) # shows contribution of each species to specific PC

  # --- combine raw data with PC scores
  df_pca <- bind_cols(df_filtered_wider_reordered,
                      as_tibble(res.pca$ind$coord), # PC1, PC2 and PC3 values
                      "Dim.1_Var" = eig.val[1,2], # % of varience explained for PC1
                      "Dim.2_Var" = eig.val[2,2],# % of varience explained for PC2
                      "Dim.3_Var" = eig.val[3,2]) # % of varience explained for PC3

  # --- calculate linear trend and change polarisation of PC1 if trend is positive
  # (8 out of 25 sites had an "inverted" trend --> just checking the polarization of first PC value still left 2 out of 25 inverted.. so I now use the linear trend to check for the overall polarization)
  change_Dim.1 <- lm(df_pca$Dim.1 ~ df_pca$Age_kaBP)
  df_pca$PC1_change <- as.numeric(change_Dim.1$coefficients[2])
  df_pca$Dim.1 <- if(df_pca$PC1_change[1] > 0) {df_pca$Dim.1 * -1} else {df_pca$Dim.1} # checks polarization of PC1 trend and if positive multiplies by -1

  df_pca
}

# --- apply PCA function (PCA.analysis) to data
df_pca <- lapply(id, PCA.analysis) %>%
  bind_rows %>% # combines the result to one df
  left_join(df_coreList[, c("ID", "lat", "lon")], by = "ID") # adds lat/lon data to df

# --- write function to run interpolate indivudal PC1 axes
PCA.interpolation <- function(id) {
  # --- filter data by ID
  df_pca <- filter(df_pca, ID == id)

  # --- interpolation of PC1 records to same time steps as set in seq.interpol (500 year bin)
  interpol.PC1 <- approx(df_pca$Age_kaBP, df_pca$Dim.1, xout = seq.interpol, rule = 1)
  interpol.PC1$y
}

# --- apply interpolation function (PCA.interpolation) to data
interpolation <- lapply(id, PCA.interpolation) %>%
  bind_cols()
colnames(interpolation) <- paste0(site, "_PC1")
interpolation <- bind_cols(t.PC1, interpolation) %>%
  dplyr::select(Age_kaBP, everything()) %>%
  filter(Age_kaBP >= 2.5 & Age_kaBP <= 23) # restrict interpolation to 2.5-23 ka (in this intervall all sites have data --> no edge issues)

# --- produce long format table out of interpolation df (needed for LOESS fit!)
interpolation_long <- interpolation %>%
  pivot_longer(-Age_kaBP, names_to = "site", values_to = "PC1_score")

# ------------------------------------------------------------------------------------------------------------------------
# --- calculate LOESS fit to PC1 scores of indivudal sites
# ------------------------------------------------------------------------------------------------------------------------
### --- Create model that will calculate LOESS fit (will do the same as under the hood in ggplot2)
model <- loess(PC1_score ~ Age_kaBP, data = interpolation_long, degree = 1, span = 0.1)
# Add predicted values from model to original dataset using broom library
interpolation_long <- augment(model, interpolation_long)

# add LOESS values to interpolation df
interpolation <- interpolation %>%
  left_join(interpolation_long[, c("Age_kaBP", ".fitted")], by = "Age_kaBP") %>%
  distinct() %>%
  rename("LOESS_fit" = ".fitted")


# ------------------------------------------------------------------------------------------------------------------------
# --- plot individual PC1 scores and LOESS fit (Figure 2a)
# ------------------------------------------------------------------------------------------------------------------------
# data frame that contains information for background shadings
df_shading <- data.frame(xmin = c(0, 11.7, 17),
                         xmax = c(11.7, 17, 24),
                         ymin = -Inf,
                         ymax = Inf,
                         name = c("Holocene", "Gloal-scale transition", "Last cold stage"),
                         col.fill = c("cadetblue4", NA, "cadetblue4")) %>%
  mutate(median_x = xmin + (xmax-xmin)/2)

plot.PC <- ggplot() +
  geom_rect(data = df_shading, # produces background shadings
            mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = df_shading$col.fill, alpha = 0.1) +
  geom_text(data = df_shading,# produces labels for shadings
            mapping = aes(x = median_x, y = ymax, label = name),
            vjust = 1.5, size = 6) + # makes the labels not "clinge" to the upper margin
  # --- interpolated PC1 scores (grey colour)
  geom_line(data = interpolation_long,
            mapping = aes(x = Age_kaBP, y = PC1_score, group = site),
            size = 0.5, colour = "grey80", linetype = 1) +
  # --- loess smooth interpolated PC1s (from ggplot2)
   geom_smooth(data = interpolation_long, mapping = aes(x = Age_kaBP, y = PC1_score),
               size = 0.5, colour = "black", linetype = 1,  method = "loess", span = 0.1, se = TRUE, method.args = list(degree=1)) +
  # # --- loess smooth interpolated PC1s (from model)
  # geom_line(data = interpolation,
  #           mapping = aes(x = Age_kaBP, y = LOESS_fit),
  #           size = 0.5, colour = "blue", linetype = 1) +
  scale_x_continuous(breaks = seq.plots, limits = c(min_age, max_age), expand = c(0.005,0.005)) +
  theme_classic(base_size = 12,
                base_line_size = 0.5) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12)) +
  labs(title = paste0("PF, interpolated PC1s (500 yr bins); LOESS fit; 2.5-23 ka"),
       x = "Age [ka BP]", y = "PC1 scores")
plot.PC

ggsave(filename = "Figure2a_PC1_LOESSfit.png", plot = plot.PC, units = "cm", height = 12, width = 30, dpi = 300)

# ------------------------------------------------------------------------------------------------------------------------
# --- plot Variance explained by PC1 on map (Figure 2b)
# ------------------------------------------------------------------------------------------------------------------------
#### --- Maps
map <- map_data("world")

# - plot
plot.var <- set.seed(1) %>% ggplot() +
  geom_polygon(data = map, aes(long, lat, group = group), fill = "grey95", colour = "grey95") +
  coord_fixed(xlim = c(-70, 25), ylim = c(-10, 75)) +
  labs(x = "Longitude [DecDeg]", y = "Latitude [DecDeg]") +
  theme(panel.background = element_rect(fill = "grey80"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "right",
        legend.direction ="vertical",
        legend.background = element_blank(),
        #legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 10),
        axis.text.y = element_text(colour = "black", size = 10),
        #plot.title = element_text(size = 14),
        legend.text = element_text(size = 10),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  geom_point(data = unique(dplyr::select(df_pca, "lon", "lat", "Dim.1_Var")), aes(x = lon, y = lat, size = Dim.1_Var),
             shape = 1) +
  scale_size_binned("PC1 variance \nexplained [%]", range = c(1, 8))
plot.var

ggsave(filename = "Figure2b_map_PC1_varianceexplained.png", plot = plot.var, units = "cm", height = 12, width = 12, dpi = 300)

# ------------------------------------------------------------------------------------------------------------------------
# --- plot global mean surface temperature (Figure 2c)
# ------------------------------------------------------------------------------------------------------------------------
# download LGMR_GMST_climo.nc from https://www.ncdc.noaa.gov/paleo/study/33112 (Osman et al., 2021)
# load temp data
Osman <- nc_open('LGMR_GMST_climo.nc')
print(Osman)
age <- ncvar_get(Osman, "age")/1000 # get age and recalculate to ka
tmp <- ncvar_get(Osman, "gmst") # get temperature
nc_close(Osman)
df_env <- tibble("Osman_Age_kaBP" = age, "Osman_global" = tmp) # save to dataframe
df_env$Osman_global <- df_env$Osman_global - mean(df_env$Osman_global[1:10]) # recalculate temperature tu mean of past 2 millennia

# --- plot SST
plot.env <- ggplot() +
  geom_rect(data = df_shading, # produces background shadings
            mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = df_shading$col.fill, alpha = 0.1) +
  geom_text(data = df_shading,# produces labels for shadings
            mapping = aes(x = median_x, y = ymax, label = name),
            vjust = 1.5, size = 6) + # makes the labels not "clinge" to the upper margin
  # --- Global mean surface temperature by Osman et al., 2021
  geom_line(data = filter(df_env, Osman_Age_kaBP >= 2.5 & Osman_Age_kaBP <= 23),
            mapping = aes(x = Osman_Age_kaBP, y = Osman_global, colour = "darkred"),
            size = 0.5, linetype = 1) +
  scale_x_continuous(breaks = seq.plots, limits = c(min_age, max_age), expand = c(0.005,0.005)) +
  theme_classic(base_size = 12,
                base_line_size = 0.5) +
  theme(legend.position = "none",
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12)) +
  #guides(col = guide_legend(nrow = 3)) +
  labs(title = paste0("Global mean surface temperature by Osman et al., 2021"),
    x = "Age [ka BP]", y = "SST anomaly [°C]")
plot.env

ggsave(filename = "Figure2c_GlobalMeanSST_Osman_etal2021.png", plot = plot.env, units = "cm", height = 12, width = 30, dpi = 300)

## ------------------------------------------------------------
# --- plot comparison of overall compositional change (LOESS fit) and global warming (temperature anomaly) --> Figure 2d
# ------------------------------------------------------------
# --- interpolation of environmental parameters at same time steps as set in seq.interpol
interpol.SST.Osman <- approx(df_env$Osman_Age_kaBP, df_env$Osman_global, xout = seq.interpol, rule = 1)

# --- combining results of EOF analysis and SST interpolations
interpolation.env <- bind_cols(t.PC1,
                               "SST_Osman" = interpol.SST.Osman$y) %>%
  filter(Age_kaBP >= 2.5 & Age_kaBP <= 23)
interpolation.env <- bind_cols(interpolation, dplyr::select(interpolation.env, -Age_kaBP))

plot.LOESSvsTemp <- ggplot(data = interpolation.env) +
  # --- LOESS vs Osman SST
  geom_point(aes(x = SST_Osman, y = LOESS_fit, colour = Age_kaBP), size = 2, shape = 16) +
  scale_color_viridis(discrete = FALSE, option = "D", direction = -1, name = "Age [ka]") +
  theme_classic(base_size = 12,
                base_line_size = 0.5) +
  theme(legend.position = c(0.1, 0.8),
        legend.title = element_text(colour = "black", size = 12),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12)) +
  labs(#title = "200 yr kernel smooth",
       x = "SST anomaly [°C]", y = "LOESS fit")
plot.LOESSvsTemp

ggsave(filename = "Figure2d_Comparison_LOESSfit_Temperature.png", plot = plot.LOESSvsTemp, units = "cm", height = 12, width = 12, dpi = 300)


# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
##### - Overview map sites (Figure 1a)
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
# download WOA18 data (statistical mean, 1/4°, annual) from https://www.ncei.noaa.gov/access/world-ocean-atlas-2018/bin/woa18.pl (Locarnini et al., 2019)
### --- read WOA18 SST data - TAKES A BIT LONG...
print("loading WOA18 data might take a while..")
SST_WOA18 <- read.csv(WOA18.data, na = na_strings, skip = 1) %>%
  clean_names() %>% # uses janitor to clean all names (all lower case with underscore)
  remove_empty(which = c("rows", "cols")) %>% # removes all empty rows and columns (sometimes produced by Excel)
  subset(select = c(1:3), latitude > -10)

### --- raster WOA18 SST data
df_raster <- SST_WOA18 %>%
  dplyr::select(x = longitude, y = latitude, z = x0)

e <- extent(c(-60, 10, -10, 80))
r <- raster(e, ncol=70, nrow=90)
r <- rasterize(df_raster[,1:2], r, df_raster[,3], fun = mean)
raster <- as.data.frame(r, xy = TRUE)

# ------------------------------------------------------------
##### plot overview map (Figure 1a)
# ------------------------------------------------------------
map <- map_data("world")

map_overview <- ggplot() +
  geom_raster(data = raster, aes(x = x, y = y, fill = layer), interpolate = TRUE) +
  geom_polygon(data = map, aes(long, lat, group = group), fill = "grey80", colour = "grey80") +
  coord_fixed(xlim = c(-60, 10), ylim = c(-10, 80), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = "grey80"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = c(0.11,0.76),
        legend.direction ="vertical",
        legend.background = element_rect(fill = "white", colour = "black"),
        axis.text.x = element_text(colour = "black", size = 10),
        axis.text.y = element_text(colour = "black", size = 10),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  geom_point(data = filter(df_coreList, plankton == plankton.group), aes(x = lon, y = lat),
             shape = 21, size = 2, colour = "black", fill = "white") +
  scale_fill_gradient2(low = "slategray3", mid = "slategray2", midpoint = 12, high = "indianred", name = "SST [°C]") # use this option for raw data maps
map_overview

ggsave(filename = "Figure1a_OverviewMap_Sites.png", plot = map_overview, units = "cm", height = 12, width = 12, dpi = 300)


# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
# --- PCA - RGB analysis (shown in Figure 1b)
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
# --- Dissimilarity metric
metric = "morisitahorn"

# --- Calculating dissimilarities for PCA
if(metric == "braycurtis")   { dcomm <- as_tibble(as(vegdist(x = df_filtered_wider_reordered[-(1:5)], method = "bray", binary = FALSE), "matrix")) }
if(metric == "morisitahorn") { dcomm <- as_tibble(as(vegdist(x = df_filtered_wider_reordered[-(1:5)], method = "horn", binary = FALSE), "matrix")) }
if(metric == "jaccard_abund"){ dcomm <- as_tibble(as(vegdist(x = df_filtered_wider_reordered[-(1:5)], method = "jaccard", binary = FALSE), "matrix")) }

# --- PCA of dissimilarities - TAKES A BIT LONG...
res.pca <- prcomp(dcomm) # runs PCA on dissimilarity matrix

# --- Getting PC projections - TAKES A BIT LONG...
pcaind <- get_pca_ind(res.pca) # package factorextra

scree.plot <- fviz_eig(res.pca, addlabels = TRUE)

# --- Copying PC projections back to the raw data
# (this is only possible because the vegdist and pca functions keep the same row order than the raw data)
df_filtered_wider_reordered$pc1 <- pcaind$coord[,1]
df_filtered_wider_reordered$pc2 <- pcaind$coord[,2]
df_filtered_wider_reordered$pc3 <- pcaind$coord[,3]

# --- Getting RGB palette for PC1,2,3
df_filtered_wider_reordered$rgb <- rgb(blue = (df_filtered_wider_reordered$pc1-min(df_filtered_wider_reordered$pc1))/(max(df_filtered_wider_reordered$pc1)-min(df_filtered_wider_reordered$pc1))*255,
                                       red = (df_filtered_wider_reordered$pc2-min(df_filtered_wider_reordered$pc2))/(max(df_filtered_wider_reordered$pc2)-min(df_filtered_wider_reordered$pc2))*255,
                                       green = (df_filtered_wider_reordered$pc3-min(df_filtered_wider_reordered$pc3))/(max(df_filtered_wider_reordered$pc3)-min(df_filtered_wider_reordered$pc3))*255,
                                       maxColorValue = 255)

# --- add lat, lon from core list
df_filtered_wider_reordered <- df_filtered_wider_reordered %>%
  left_join(df_coreList[, c("ID", "lat", "lon")], by = "ID")



# ------------------------------------------------------------------------------------------------------------------------
# --- plot RGB values in Hovmoller-like plot (Figure 1b)
# ------------------------------------------------------------------------------------------------------------------------
df_filtered_wider_reordered$lat_round <- 2.5*round(df_filtered_wider_reordered$lat/2.5) # you have to play around with the number here to see which rounding works best for you
df_filtered_wider_reordered$age_round <- 1*round(df_filtered_wider_reordered$Age_kaBP/1) # same as above, for example, substitude 1 for 3 to see plot

df_filtered_wider_reordered$pc1_norm <- (df_filtered_wider_reordered$pc1-min(df_filtered_wider_reordered$pc1))/(max(df_filtered_wider_reordered$pc1)-min(df_filtered_wider_reordered$pc1))
df_filtered_wider_reordered$pc2_norm <- (df_filtered_wider_reordered$pc2-min(df_filtered_wider_reordered$pc2))/(max(df_filtered_wider_reordered$pc2)-min(df_filtered_wider_reordered$pc2))
df_filtered_wider_reordered$pc3_norm <- (df_filtered_wider_reordered$pc3-min(df_filtered_wider_reordered$pc3))/(max(df_filtered_wider_reordered$pc3)-min(df_filtered_wider_reordered$pc3))

# --- group by 1) core + age, 2) lat + age
# (PC scores are averaged according 1) and 2) so for later tile plotting there is one RGB colour per tile, otherwise all RGB colour tiles are stacked on top of each other)
df_filtered_group <- df_filtered_wider_reordered %>%
  group_by(Site, age_round) %>% # grouped by site and age
  summarise(pc1_norm_mean = mean(pc1_norm),
            pc2_norm_mean = mean(pc2_norm),
            pc3_norm_mean = mean(pc3_norm)) %>%
  inner_join(df_filtered_wider_reordered[, c("Site", "lat_round")], by = "Site") %>%
  distinct() %>% # remove duplicate rows
  group_by(lat_round, age_round) %>% # group a second time, but this time by lat and age
  summarise(pc1_norm_mean = mean(pc1_norm_mean),
            pc2_norm_mean = mean(pc2_norm_mean),
            pc3_norm_mean = mean(pc3_norm_mean))

# --- Grid RGB values
plot.hovmoller.gridded <-
  ggplot(data = df_filtered_wider_reordered) +
  # - RGB tiles averaged (here RGB colours are averaged by tile; one colour per tile)
  geom_tile(data = df_filtered_group,
            aes(x = age_round, y = lat_round,
                fill = rgb(blue = pc1_norm_mean, # changing the order around to change color scheme
                           red = pc2_norm_mean,
                           green = pc3_norm_mean))) +
  scale_fill_identity() +
  geom_point(aes(x = Age_kaBP, y = lat), colour = "grey80", size = 0.5) +
  scale_colour_identity() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "none",
        legend.direction ="vertical",
        legend.background = element_blank(),
        legend.title = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 9),
        axis.title = element_text(colour = "black", size = 12),
        axis.text = element_text(colour = "black", size = 10),
        plot.title = element_text(size = 12)) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_x_continuous(breaks = seq(0,24, by = 3), expand=c(0,0)) +
  labs(x = "Age [ka]", y = "Latitude")
plot.hovmoller.gridded

ggsave(filename = "Figure1b_AssemblageComposition.png", plot = plot.hovmoller.gridded, units = "cm", height = 12, width = 16, dpi = 300)

##### ---- colorblind-friendly version of Figure 1b -----
# --- color-code each composition with mix() to generate a colorblind-friendly version of Figure 1b
library(colormod)

df_filtered_group$amount_pc1 <- df_filtered_group$pc1_norm_mean/(df_filtered_group$pc1_norm_mean + df_filtered_group$pc2_norm_mean + df_filtered_group$pc3_norm_mean)
df_filtered_group$amount_pc2 <- df_filtered_group$pc2_norm_mean/(df_filtered_group$pc1_norm_mean + df_filtered_group$pc2_norm_mean + df_filtered_group$pc3_norm_mean)
df_filtered_group$amount_pc3 <- df_filtered_group$pc3_norm_mean/(df_filtered_group$pc1_norm_mean + df_filtered_group$pc2_norm_mean + df_filtered_group$pc3_norm_mean)

df_filtered_group <- df_filtered_group %>%
  rowwise() %>%
  mutate(
    # # individual colours chosen from viridis palette (https://waldyrious.net/viridis-palette-generator/)
    # rgb_pc1 = mix('#fde725', '#31688e', pc1_norm_mean), # colour 1: yellow, colour 2: blue
    # rgb_pc2 = mix('#35b779', '#440154', pc2_norm_mean), # colour 1: green , colour 2: purple
    # individual colours chosen from plasma palette (https://waldyrious.net/viridis-palette-generator/)
    rgb_pc1 = mix('#f0f921', '#9c179e', pc1_norm_mean), # colour 1: yellow, colour 2: purple
    rgb_pc2 = mix('#ed7953', '#0d0887', pc2_norm_mean), # colour 1: orange, colour 2: dark blue
    rgb_mix = mix(rgb_pc1, rgb_pc2, pc3_norm_mean)
  )

# --- Grid RGB values
plot.hovmoller.gridded.col <-
  ggplot(data = df_filtered_wider_reordered) +
# - RGB tiles via mix() and shade()/tint() from colormod
geom_tile(data = df_filtered_group,
          aes(x = age_round, y = lat_round,
              fill = rgb_mix)) +
  scale_fill_identity() +
  geom_point(aes(x = Age_kaBP, y = lat), colour = "grey80", size = 0.5) +
  scale_colour_identity() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "none",
        legend.direction ="vertical",
        legend.background = element_blank(),
        legend.title = element_text(colour = "black", size = 10),
        legend.text = element_text(size = 9),
        axis.title = element_text(colour = "black", size = 12),
        axis.text = element_text(colour = "black", size = 10),
        plot.title = element_text(size = 12)) +
  scale_y_continuous(expand=c(0, 0)) +
  scale_x_continuous(breaks = seq(0,24, by = 3), expand=c(0,0)) +
  labs(x = "Age [ka]", y = "Latitude"
  )
plot.hovmoller.gridded.col

ggsave(filename = "Figure1b_AssemblageComposition_colourblind_friendly.png", plot = plot.hovmoller.gridded.col, units = "cm", height = 12, width = 16, dpi = 300)
