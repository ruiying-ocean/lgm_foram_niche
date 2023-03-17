library(tidyverse)
library(vegan)
foram_abun <- read_csv("data/Strack_2022/deglacial_foram_with_sst.csv")
foram_abun <- foram_abun %>% mutate(agebin = round(Age_kaBP)) %>% select(c(Site, Lat, Lon, Age_kaBP, agebin, SST, Species, Rel_abundance))
foram_abun <- foram_abun %>% filter(!is.na(Rel_abundance))

# to wide table
foram_abun_wide <- foram_abun %>%
  pivot_wider(names_from = Species, values_from = Rel_abundance, values_fill = 0)

# group species
# identify different names
# species <- foram_abun$Species %>% unique()
# in_the_table <- species %in% symbiosis_tbl$Species
# species[!in_the_table]

foram_abun <- foram_abun %>% mutate(Species=
                                      recode(Species,
                                             "Globigerinoides_ruber_ruber" =  "Globigerinoides_ruber",
                                             "Globigerinoides_ruber_albus" = "Globigerinoides_ruber",
                                             "Globigerinoides_ruber_ruber_and_Globigerinoides_ruber_albus" = "Globigerinoides_ruber",
                                             "Globorotalia_menardii_and_Globorotalia_tumida" = "Globorotalia_menardii"))

source("code/data_clean/read_symbiosis_table.R")
symbiosis_tbl$Species <- gsub(" ", "_",  symbiosis_tbl$Species)
symbiosis_tbl <- symbiosis_tbl %>% select(!short_name)
foram_abun <- left_join(foram_abun, symbiosis_tbl)
grouped_foram_abun <- foram_abun %>% group_by(Site, Lat, Lon, Age_kaBP, Symbiosis) %>% 
  reframe(Rel_abundance = sum(Rel_abundance), SST=SST) %>% filter(Symbiosis != "Undetermined")
grouped_foram_abun <- grouped_foram_abun %>% distinct()

grouped_foram_abun <- grouped_foram_abun %>% 
  mutate(lat_group = case_when(
    Lat < 50 & Lat > 35 ~ "temperate",
    Lat < 35 ~ "tropical",
    Lat > 50 ~ "(sub)polar"
  ))

# optimal temperature is not able to calculate because it is relative abundance
# opt_foram_abun <- foram_abun %>% group_by(agebin, Species) %>% reframe(Rel_abundance = max(Rel_abundance))
# opt_foram_abun <- opt_foram_abun %>% filter(Rel_abundance > 0)
# opt_foram_temp <- left_join(opt_foram_abun, foram_abun, by=c("agebin", "Species", "Rel_abundance"))
# opt_foram_temp %>% ggplot(aes(x=agebin, y=SST)) + geom_point() + facet_wrap(~Species)

site_list <- foram_abun_wide$Site %>% unique()
n = 1603 # number from test
site <- array(numeric(),n)
mid_age <- array(numeric(),n)
sst1 <- array(numeric(),n)
sst2 <- array(numeric(),n)
lat <- array(numeric(),n)
delta_sst <- array(numeric(),n)
delta_dis <- array(numeric(),n)
delta_rel <- array(numeric(),n)
accum_dis <- array(numeric(),n)
accum_sst <- array(numeric(),n)
accum_rel <- array(numeric(),n)

## start loop
count=0
for (i in site_list) {
  # subset each site
  site_data <- foram_abun_wide %>% filter(Site==i) %>% arrange(Age_kaBP)
  
  site_id <- site_data %>% pull(Site) %>% unique()
  site_lat <- site_data %>% pull(Lat) %>% unique()
  site_age <- site_data %>% pull(Age_kaBP)
  site_sst <- site_data %>% pull(SST)  
  site_age_length <- length(site_age)
  site_abun <-  grouped_foram_abun %>% 
    filter(Site==site_id, Symbiosis=="Yes") %>%
    pull(Rel_abundance)
  site_species <- site_data %>% select(Neogloboquadrina_dutertrei:Turborotalita_humilis)
  
  for(j in seq(site_age_length)){
    count = count + 1
    max_age = max(site_age_length)

    if (j == max_age) {
      delta_dis[count] = NA
      delta_sst[count] = NA
      mid_age[count] = NA
      sst1[count] = NA
      sst2[count] = NA
      accum_dis[count]= NA
      lat[count] = NA
      accum_sst[count]= NA
      delta_rel[count] = NA
      accum_rel[count] = NA
    } else {
      
    horn <- rbind(site_species[j, ], site_species[j+1,]) %>% 
      vegdist(., method="horn", na.rm = FALSE)
    
    accum_horn <- rbind(site_species[j, ], site_species[max_age,]) %>% 
      vegdist(., method="horn", na.rm = FALSE)
    
    dtmp <- site_sst[j] - site_sst[j+1]
    accum_temp <- site_sst[j] - site_sst[max_age]
    
    # saving data
    mid_age[count] = (site_age[j]+site_age[j+1])/2
    sst1[count]  <- site_sst[j]
    sst2[count]  <- site_sst[j+1]
    
    delta_sst[count] = dtmp
    delta_dis[count] = horn
    accum_sst[count] = accum_temp
    accum_dis[count] = accum_horn
    lat[count] = site_lat
    site[count] = site_id
    
    # delta relative abundance
    delta_rel[count] = site_abun[j]-site_abun[j+1]
    
    # delta relative abundance
    accum_rel[count] = site_abun[j]- site_abun[max_age]
    }
  }
}

## plot variable: assemblage composition index; relative abundance; 
## plot type: (1) nearest; (2) accumulative
## plot option: (1) time series; (2) vs sst

df <- data.frame(site=site, mid_age = mid_age, 
                 sst1=sst1, sst2=sst2, lat=lat,
                 dsst = delta_sst, dsim = delta_dis, dabun = delta_rel,
                 accum_sst = accum_sst, accum_dis = accum_dis, accum_abun = accum_rel)
df <- df %>%
  mutate(lat_group = case_when(
    lat < 50 & lat > 35 ~ "temperate",
    lat < 35 ~ "tropical",
    lat > 50 ~ "(sub)polar"
  ))
# First, local temperature are differently influenced
p1 <- foram_abun %>% filter(!is.na(SST)) %>% group_by(group= Site) %>% mutate(min_SST = min(SST, na.rm=T)) %>% 
  mutate(SST_anomaly = SST- min_SST) %>%
  ggplot(aes(x=Age_kaBP, y =SST_anomaly)) + geom_line(linewidth=1, aes(lty=Site, color=Lat))+
  theme_linedraw() + scale_color_viridis_c(direction=-1) + guides(lty="none") + 
  theme(legend.position = c(.85, .7), legend.background = element_blank())

## second, aseemblage change is mostly happening in transitional areas
df %>% ggplot(aes(x=mid_age, y=accum_abun)) + geom_line(aes(color=as.character(lat)), linewidth=1) + facet_wrap(~lat)
df %>% ggplot(aes(x=mid_age, y=accum_dis)) + geom_line(aes(color=as.character(lat)), linewidth=1) + facet_wrap(~lat)

# Third, this change is not influenced by SST
df %>% ggplot(aes(x=delta_sst, y=delta_dis)) +  geom_point(aes(color=lat_group))
df %>% ggplot(aes(x=delta_sst, y=delta_rel)) +  geom_point(aes(color=lat_group))

## conclusion
## (1) global warming is impacting ocean differently, the polar ocean will experience highest anamoly
## (2) foram in polar/tropical will stay (little compositional change), transitional area will see most dispersal
## (3) the change is not driven by temperature

# library(patchwork)
# p <- p1+p2+p3
# p
# ggsave("output/degalcial_composition.pdf", width=10, height=4)
