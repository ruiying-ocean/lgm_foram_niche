library(tidyverse)
library(patchwork)
library(hrbrthemes)
library(RColorBrewer)
library(ncdf4)
library(zoo)

## --------------- Observational data -------------------
# temperature and CO2 data source: https://ourworldindata.org/co2-and-other-greenhouse-gas-emissions
past_temp <- read_csv("data/model_drived//observed_global_temperature_anamoly.csv")
past_co2 <- read_csv("data/model_drived/observed_co2.csv")
names(past_temp)[2] <- "Global_Air_Temp_Anomaly"

### AMOC observations, source: https://www.bodc.ac.uk/data/published_data_library/catalogue/10.5285/e91b10af-6f0a-7fa7-e053-6c86abc05a09
# nc_data <- nc_open('data/modern_env/RAPID/moc_transports_200404_20201215.nc')
# amoc_transport <- ncvar_get(nc_data, "moc_mar_hc10")
# rapid_time <-  ncvar_get(nc_data, "time") #time in days since 2004-04-01
# rapid_inittime <- "2004-04-01" %>% as.Date(., format = "%Y-%m-%d")
# rapid_time <- rapid_inittime + rapid_time
# nc_close(nc_data) 

### Global dust loading, from Kok et al. 2023 Nature Review
dust_nc <- nc_open("data/modern_env/DustCOMM_global_historical_evolution_dust_loading.nc")
gbl_dust_loading <- ncvar_get(dust_nc, "mean")
dust_loading_time <- ncvar_get(dust_nc, "end_year")
nc_close(dust_nc) 
gbl_dust_loading <- tibble(year = dust_loading_time, global_dust_loading = gbl_dust_loading) %>%
  mutate(year = as.Date(as.yearmon(year)))

### construct AMOC dataframe
# rapid_amoc <- tibble(date = rapid_time, amoc_transport = amoc_transport)
# rapid_amoc <- rapid_amoc %>% mutate(rapid_amoc_yearly = rollmeanr(rapid_amoc$amoc_transport, 30*12, fill = NA))

## Longer term AMOC, but proxy
## source: https://www.nature.com/articles/s41586-018-0007-4#MOESM5
amoc_proxy <- read_csv("data/modern_env/Long_term_AMOC_proxy.csv") %>%
  mutate(transport =`Tsub AMOC dipole`*2.3,
         Year = as.Date(as.yearmon(Year)))

## --------------- Model data -------------------
model_pco2_temp <- read_csv("data/model_drived/model_pCO2_temperature.csv")
#model_pco2_temp <- model_pco2_temp %>% mutate(Age_group = as.factor(if_else(Year < 2022, "History", "Future")))
model_amoc <- read_csv("data/model_drived/model_moc.csv") %>% mutate(Year = as.Date(as.yearmon(Year)))

model_pco2_temp$Year <- model_pco2_temp$Year - 0.5
model_pco2_temp$Group <- as.factor(model_pco2_temp$Group)

# as the data, set average 1960-1990 as standard
Average_standard <- model_pco2_temp %>% filter(Year > 1960 & Year < 1990, Group=='PI') %>% summarise(avg_temp = mean(Temperature))
model_pco2_temp <- model_pco2_temp %>% group_by(Group) %>% mutate(delta_T= Temperature - Average_standard$avg_temp) %>% ungroup()

## AMOC baseline for both model and data (2004 as baseline)
amoc_proxy_baseline <- amoc_proxy %>% filter(Year==as.Date("1775-01-01", format = "%Y-%m-%d")) %>% pull(transport)
model_amoc_baseline <- model_amoc %>% filter(Year==as.Date("1765-07-01", format = "%Y-%m-%d")) %>% pull(AMOC_max)

## --------------- Plot data -------------------
p_temp <- ggplot() + geom_line(data=model_pco2_temp, aes(x=Year, y=delta_T, group=Group, color=Group), linewidth=1.2) + 
  geom_line(data=past_temp, aes(x=Year, y=Global_Air_Temp_Anomaly)) +
  theme_ipsum() + scale_color_ipsum() + 
  labs(color="Global warming by 2100\nrelative to 1760 (°C)", x='', y="Tempearature Relative to 1960-1990 (°C)") + 
  theme(legend.position=c(0.3, 0.78))

p_co2 <- ggplot() + geom_line(data=model_pco2_temp, aes(x=Year, y=`pCO2 (ppm)`, group=Group, color=Group), linewidth=1.2) + 
  geom_point(data=past_co2, aes(x=Year, y=CO2_ppm), size=1.5) +
  theme_ipsum() + scale_color_ipsum()  + labs(color="Global warming by 2100\nrelative to 1760 (°C)", x='', y="pCO2 (ppm)") +
  theme(legend.position=c(0.3, 0.78))

p_amoc <- ggplot() + geom_line(data=model_amoc, aes(x=Year, y=AMOC_max-model_amoc_baseline, group=Group, color=Group), linewidth=1.2)+
  theme_bw()+
  theme(legend.position="none") + geom_line(data=amoc_proxy, aes(x=Year, y=transport-amoc_proxy_baseline)) +
  scale_x_date(limits = as.Date(c("1765-01-01","2100-01-01")))+
  labs(y="Sv", title="Modern AMOC transport", x="") + 
  scale_color_manual(values=c("#FEEDDE", "#FDBE85", "#FD8D3C", "#D94701", "#66C2A5")) +
  theme(legend.position=c(0.4, 0.85), legend.direction = "horizontal", legend.background = element_blank(), legend.title = element_blank())+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(text=element_text(family="Fira Sans", size=font_size))

fig_futureco2 <- ggplot() + geom_line(data=model_pco2_temp, aes(x=Year, y=`pCO2 (ppm)`, group=Group, color=Group), linewidth=1.2) + 
  geom_point(data=past_co2, aes(x=Year, y=CO2_ppm), size=1.5) + theme_bw() +
  labs(y="ppm", title=TeX("Modern $pCO_2$")) + scale_color_manual(values=c("#FEEDDE", "#FDBE85", "#FD8D3C", "#D94701", "#66C2A5"))+
  theme(legend.position=c(0.4, 0.85), legend.direction = "horizontal", legend.background = element_blank(), legend.title = element_blank())+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  theme(text=element_text(family="Fira Sans", size=font_size))

fig_dust_loading <- ggplot() + geom_line(data=gbl_dust_loading, aes(x=year, y=global_dust_loading)) + 
  theme_bw() + xlim(c(1765, 2100))+
  labs(y="Tg", title="Modern global dust loading", x="Year")+
  scale_x_date(limits = as.Date(c("1765-01-01","2100-01-01")))+
  theme(text=element_text(family="Fira Sans", size=font_size))

p <- p_temp +  p_co2 + plot_annotation(tag_levels = 'A')
# png("./output/temperature_gradient.jpg", res=300, width = 10, height = 5, unit="in")
# print(p)
# dev.off()
