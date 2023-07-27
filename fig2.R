library(ggpubr)

antell_2021 <- readRDS("data/Antell_realaysis.RDS")
symbiosis_tbl <- read_csv("~/foram_core/fg/foram_sp_db.csv") %>% 
  mutate(sp = map_vec(Name, species_abbrev)) %>% mutate(fg=case_when(
    Symbiosis=="No" & Spinose == "No"~"Symbiont-barren Non-Spinose",
    Symbiosis=="Yes" & Spinose == "Yes"~"Symbiont-obligate Spinose",
    Symbiosis=="No" & Spinose == "Yes"~"Symbiont-barren Spinose",
    Symbiosis=="Yes" & Spinose == "No"~"Symbiont-facultative Spinose",
  ))
antell_2021 <- merge(antell_2021,symbiosis_tbl, by="sp")

# Load the ncdf4 library
library(ncdf4)

# Open the NetCDF file
nc <- nc_open("~/lgm_carbon/data/Tierney2020_DA_ocn_regrid.nc")

# Read the "deltaSST" variable
deltaSST <- ncvar_get(nc, "deltaSST")
lat <- ncvar_get(nc, "lat")

# Calculate mean of different latitudes
lower_lat_mean <- mean(deltaSST[abs(lat) <= 45, ], na.rm = TRUE)
higher_lat_mean <- mean(deltaSST[abs(lat) >= 45, ], na.rm = TRUE)
global_mean <- mean(deltaSST, na.rm = TRUE)

# Close the NetCDF file
nc_close(nc)

## LGM species-level
sp_opt <- thermal_opt(obs_sp_smooth)

select_species <- c("G. bulloides", "N. pachyderma",
                    "N. dutertrei", "N. incompta", "G. inflata", "G. glutinata",
                    "G. ruber albus", "G. ruber ruber","G. truncatulinoides",
                    "P. obliquiloculata", "T. quinqueloba","T. sacculifer")

sst_change <- c(global_mean, higher_lat_mean, lower_lat_mean, global_mean, global_mean,
                global_mean, lower_lat_mean, lower_lat_mean, lower_lat_mean, lower_lat_mean,
                higher_lat_mean, lower_lat_mean)

species_sst_change <- data.frame(species = select_species, delta_sst = sst_change*-1) %>%
  left_join(sp_opt, by = "species",multiple = "all") %>% group_by(species)%>%
  mutate(delta_opt = lag(model_x)- model_x)


fig2a <- ggscatter(data=antell_2021,
                   x='delta_temp', y='delta_pe', alpha=0.1, size=2,
                   ggtheme = theme_bw(), xlim=c(-5,5),
                   add="reg.line",
                   cor.coef=TRUE)

ggscatter(data=species_sst_change, color="latitude",
          x='delta_sst', y='delta_opt', alpha=1, size=2,
          ggtheme = theme_bw(), xlim=c(-5,5))

## high lat -> genetic diversity
ggscatter(data=jonkers_2019, color="y",
          x='real.dT', y='dSST.niche.evolution', alpha=1, size=2,
          ggtheme = theme_bw(), xlim=c(-5,5))+scale_colour_viridis_b()

fig2a + 
  geom_point(data=species_sst_change, aes(x = delta_sst, y = delta_opt), shape=4, color="blue")+
  geom_point(data=jonkers_2019, aes(y=dSST.niche.evolution, x=real.dT), shape=5, color="red")+
  geom_smooth(data=jonkers_2019, aes(y=dSST.niche.evolution, x=real.dT), method="lm",se=F,color="red")+
  geom_abline(slope=0)
  
fig2_model <- genie_fg_raw %>% group_by(age) %>% summarise(mean_sst = mean(sst)) %>%
  right_join(thermal_opt(genie_fg_smooth), multiple = "all")

fig2_model <- fig2_model %>%
  group_by(species) %>% 
  mutate(delta_sst = mean_sst - min(mean_sst[age == "pi"]),
         pi_opt = min(model_x[age == "pi"]),
         pi_abundance =max_y[age == "pi"]) %>%
  mutate(delta_opt = model_x - pi_opt,
         delta_abundance = max_y/pi_abundance)

fig2b <- ggplot() +  
  geom_point(data=fig2_model, aes(x = delta_sst, y = delta_opt, color = species)) +
  geom_line(data=fig2_model, aes(x = delta_sst, y = delta_opt, color = species),linewidth=1)+
  theme(legend.position = "none") +
  xlim(-5,5)+ylim(-21,21)

fig2a + fig2b
## RY re-alaysis
## two location, two times (4 SSTs)
##      LocA       LocB
## PI   A0(5 °C)   B0(4 °C)
## Now  A1(6 °C)   B1(5 °C)

## TS.SST -> pre-industrial sst of trap locaton (B0)
## TS.real.SST -> in situ collection SST (B1)
## real.dT -> in situ warming (TS.real.SST - TS.SST) (B1- B0)
## dSST.most.sim.lt -> SST.most.sim.lt (PI) - TS.SST (PI) (A0-B0)
## trap.trend, positive or negative of dSST.most.sim.lt
## real.trend, positive or negative of real.dT
## if consistent, they are moving to the "right" place

## what I want to test
## if B1 = A0, then niche stable
## if B1 > A0, then niche plasticity
## if B1 < A0, then reverse
## plot, x= (B1-B0), y=(B1-A0)
## niche evolution/temperature change

jonkers_2019 <- readRDS("data/plot_niche_Jonkers.RDS")
jonkers_2019 <-jonkers_2019 %>% filter(consistent=TRUE)


# fig3b <- fig3b+
#   annotate(geom = "polygon", x = c(-Inf, 0, 0, -Inf), y = c(-5, 0, -Inf, -Inf), fill = "red", alpha = 0.1 )+
#   annotate(geom = "polygon", x = c(0, 0, Inf, Inf), y = c(0, Inf, Inf,3.5), fill = "red", alpha = 0.1 )
# 
# fig3b <- fig3b+ annotate(geom = "polygon", x = c(-Inf, 0, -Inf), y = c(-5, 0, 0), fill = "grey", alpha = 0.1 )+
#    annotate(geom = "polygon", x = c(Inf, 0, Inf), y = c(0, 0, 3.5), fill = "grey", alpha = 0.1 )
# 
# fig3b <- fig3b+ annotate(geom = "polygon", x = c(-Inf, 0, 0,-Inf), y = c(0, 0, Inf,Inf), fill = "blue", alpha = 0.1 )+
#   annotate(geom = "polygon", x = c(Inf, 0, 0,Inf), y = c(0, 0, -Inf,-Inf), fill = "blue", alpha = 0.1 )
# 
# fig3b <- fig3b+xlim(-9,9)+ylim(-9,9)
# fig3b <- fig3b + labs(x="Temperature change relative to pre-industrial age", y="optimal temperature change")
# 
# grid_data <- expand.grid(x = 1:10, y = 1:10)
# fig3a <- ggplot(grid_data, aes(x, y))+
#   annotate(geom = "polygon", x = c(-Inf, 0, 0, -Inf), y = c(-5, 0, -Inf, -Inf), fill = "red", alpha = 0.1 )+
#   annotate(geom = "polygon", x = c(0, 0, Inf, Inf), y = c(0, Inf, Inf,3.5), fill = "red", alpha = 0.1 )+
#   annotate(geom = "polygon", x = c(-Inf, 0, -Inf), y = c(-5, 0, 0), fill = "grey", alpha = 0.1 )+
#   annotate(geom = "polygon", x = c(Inf, 0, Inf), y = c(0, 0, 3.5), fill = "grey", alpha = 0.1 )+
#   annotate(geom = "polygon", x = c(-Inf, 0, 0,-Inf), y = c(0, 0, Inf,Inf), fill = "blue", alpha = 0.1 )+
#   annotate(geom = "polygon", x = c(Inf, 0, 0,Inf), y = c(0, 0, -Inf,-Inf), fill = "blue", alpha = 0.1 )+
#   xlim(-9,9)+ylim(-9,9)
# 
# fig3a <- fig3a + annotate("text", x=-5,y=5,label="no acclimation + cold ecology", color="blue")+
#    annotate("text", x=5,y=-5,label="no acclimation + cold ecology", color="blue")
# fig3a <- fig3a + annotate("text", x=5,y=5,label="full acclimation + warm ecology", color="red")+
#  annotate("text", x=-5,y=-5,label="full acclimation + warm ecology", color="red")
# fig3a <- fig3a + annotate("text", x=6,y=1,label="imperfect\nacclimation", color="grey")+
#    annotate("text", x=-6,y=-1,label="imperfect\nacclimation", color="grey")
# fig3b=fig3b+theme(legend.position = "none")


fig2a
fig2b
fig2c
