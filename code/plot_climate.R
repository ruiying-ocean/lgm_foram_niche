library(tidyverse)
library(latex2exp)
library(patchwork)
library(ncdf4)
library(thematic)
library(showtext) # for the font
library(cmocean)

## Data sources
## --------------- Deglacial Obs -------------------
## Deglacial pCO2
deg_pco2 <- read_delim("data/deglaciation_climate/CO2_stack_156K_spline_V2.txt",quote ="\t")

## Deglacial pCO2
#isotopes <- read_delim("data/deglaciation_climate/Lippold-etal_2019_isotopes.tab")
deg_amoc <- read_csv("data/deglaciation_climate/deglacial_AMOC.csv") %>% mutate(`Age [ka BP]` = `Age (BP)`/1000)
deg_iron <- read_delim("data/deglaciation_climate/EDC_DustFlux_25yr.tab")

## --------------- Modern Obs -------------------
## Global dust loading, from Kok et al. 2023 Nature Review
dust_nc <- nc_open("data/modern_env/DustCOMM_global_historical_evolution_dust_loading.nc")
dust_sigma <- ncvar_get(dust_nc, "pos1sigman")
gbl_dust_loading <- ncvar_get(dust_nc, "mean")
dust_loading_time <- ncvar_get(dust_nc, "end_year")
nc_close(dust_nc)
remove(dust_nc)
gbl_dust_loading <- tibble(year = dust_loading_time, global_dust_loading = gbl_dust_loading, sigma=dust_sigma) %>%
    mutate(year = as.Date(as.yearmon(year)))

# historical temperature and CO2 data source: https://ourworldindata.org/co2-and-other-greenhouse-gas-emissions
hist_temp <- read_csv("data/model_drived/observed_global_temperature_anamoly.csv")
hist_co2 <- read_csv("data/model_drived/observed_co2.csv")
names(hist_temp)[2] <- "Global_Air_Temp_Anomaly"

## Late Holocene AMOC
## source: https://www.nature.com/articles/s41586-018-0007-4#MOESM5
amoc_proxy <- read_csv("data/modern_env/Long_term_AMOC_proxy.csv") %>%
  mutate(transport =`Tsub AMOC dipole`*2.3,
         Year = as.Date(as.yearmon(Year)),
         SE = `2SE`/2)

### AMOC observations, source: https://www.bodc.ac.uk/data/published_data_library/catalogue/10.5285/e91b10af-6f0a-7fa7-e053-6c86abc05a09
# nc_data <- nc_open('data/modern_env/RAPID/moc_transports_200404_20201215.nc')
# amoc_transport <- ncvar_get(nc_data, "moc_mar_hc10")
# rapid_time <-  ncvar_get(nc_data, "time") #time in days since 2004-04-01
# rapid_inittime <- "2004-04-01" %>% as.Date(., format = "%Y-%m-%d")
# rapid_time <- rapid_inittime + rapid_time
# nc_close(nc_data) 

# rapid_amoc <- tibble(date = rapid_time, amoc_transport = amoc_transport)
# rapid_amoc <- rapid_amoc %>% mutate(rapid_amoc_yearly = rollmeanr(rapid_amoc$amoc_transport, 30*12, fill = NA))

## --------------- Model data -------------------
model_pco2_temp <- read_csv("data/model_drived/model_pCO2_temperature.csv")
#model_pco2_temp <- model_pco2_temp %>% mutate(Age_group = as.factor(if_else(Year < 2022, "History", "Future")))
model_amoc <- read_csv("data/model_drived/model_moc.csv") %>% mutate(Year = as.Date(as.yearmon(Year)))

model_pco2_temp$Year <- model_pco2_temp$Year - 0.5
model_pco2_temp$Group <- factor(model_pco2_temp$Group, levels=c('PI', '1p5', '2', '3' ,'4'))
model_amoc$Group <- factor(model_amoc$Group, levels=c('PI', '1p5', '2', '3' ,'4'))

# as the data, set average 1960-1990 as standard
Average_standard <- model_pco2_temp %>% filter(Year > 1960 & Year < 1990, Group=='PI') %>% summarise(avg_temp = mean(Temperature))
model_pco2_temp <- model_pco2_temp %>% group_by(Group) %>% mutate(delta_T= Temperature - Average_standard$avg_temp) %>% ungroup()

## AMOC baseline for both model and data (2004 as baseline)
amoc_proxy_baseline <- amoc_proxy %>% filter(Year==as.Date("1750-01-01", format = "%Y-%m-%d")) %>% pull(transport)
model_amoc_baseline <- model_amoc %>% filter(Year==as.Date("1765-07-01", format = "%Y-%m-%d")) %>% pull(AMOC_max)

## --------------- Plot deglacial warming -------------------
theme_set(theme_classic())
thematic_on(font = font_spec("Fira Sans", scale=1))
pal <- c("#00798c", "#d1495b", "#edae49","#F07818","#6495ED", "#A155B9")
model_pal <- c("gray70","#FEEDB0FF","#F29567FF","#CE4356FF","#831C63FF")

add_chronology <- function(p, add_text=FALSE, label_position=400){
  p <- p +
    annotate("rect", xmin =14.7, xmax = 18.1, ymin = -Inf, ymax = Inf, fill = "grey", alpha = .15, color = NA)+
    annotate("rect", xmin =12.9, xmax = 11.6, ymin = -Inf, ymax = Inf, fill = "grey", alpha = .15, color = NA)+
    annotate("rect", xmin = 0.2, xmax = -0.2, ymin = -Inf, ymax = Inf, fill = "grey", alpha = .15, color = NA)+
    
  if (add_text){
    p <- p +annotate("text", x =13.8, y=label_position, label="BA")+
      annotate("text", x =12.2, y=label_position, label="YD")+
      annotate("text", x =19.5, y=label_position, label="LGM")+
      annotate("text", x =16.3, y=label_position, label="HS1")+
      annotate("text", x =5, y=label_position, label="Holocene")
  }
  return(p)
}

# global co2
fig_co2 <- deg_pco2 %>% filter(`Age [ka BP]` <=20) %>%  
  ggplot(., aes(x=`Age [ka BP]`)) +
  geom_line(aes(y=`CO2 [µmol/mol]`), color=pal[1]) +
  geom_ribbon(aes(ymin=`CO2 [µmol/mol]` - `CO2 std dev [±]`,
                  ymax=`CO2 [µmol/mol]` + `CO2 std dev [±]`),
                  alpha=.2, fill=pal[1]) +
  labs(y="ppm", x="") +
  xlim(c(20, -0.2)) + 
  theme(legend.position = "none")

fig_co2 <- fig_co2 + annotate("rect", xmin = 0.5, xmax=-0.2, ymin=280, ymax=410,
                fill="white", color="black", linetype="dotted", alpha=0)

# Pa/Th, lower value -> vigorous NADW 
fig_amoc <- deg_amoc %>% filter(`Age [ka BP]` <=20) %>%
  ggplot(., aes(x=`Age [ka BP]`)) +
  geom_point(aes(y=`Individual Pa/Th data points`), size=.5, alpha=0.4) +
  geom_ribbon(aes(ymin=`Pa/Th composite curve (9-point moving average)` - `2 standard error`/2,
                  ymax=`Pa/Th composite curve (9-point moving average)` + `2 standard error`/2),
              alpha=.2, fill=pal[2])+
  geom_line(aes(y=`Pa/Th composite curve (9-point moving average)`), color=pal[2])+
     xlim(c(20, -0.2))+ylim(c(0.09, 0.045))+
  labs(y=TeX("$^{231}$Pa/$^{230}Th$"),x="")

## Add annotation
fig_amoc <- fig_amoc + geom_segment(aes(x = 2, y = 0.08, xend = 2, yend = 0.07),
                      arrow = arrow(length = unit(0.06, "inches"), ends = "both")) +
  annotate("text", x = 2, y = 0.068, label = "Strong")+
  annotate("text", x = 2, y = 0.082, label = "Weak")

# Dome C core deg_iron
fig_iron <- deg_iron %>% filter(`Age [ka BP]`<=20) %>% ggplot(., aes(x=`Age [ka BP]`)) +
  geom_line(aes(y=`Dust flux [mg/m**2/a]`), color=pal[3]) + xlim(c(20, -0.2)) + 
  labs(y=TeX("mg m$^{-2}$ yr$^{-1}$"), x="Thousands of Years (BCE)")

## --------------- Plot modern/future warming -------------------

# p_temp <- ggplot() + geom_line(data=model_pco2_temp, aes(x=Year, y=delta_T, group=Group, color=Group), linewidth=1.2) + 
#   geom_line(data=hist_temp, aes(x=Year, y=Global_Air_Temp_Anomaly)) +
#   theme_ipsum() + scale_color_ipsum() + 
#   labs(color="Global warming by 2100\nrelative to 1760 (°C)", x='', y="Tempearature Relative to 1960-1990 (°C)") + 
#   theme(legend.position=c(0.3, 0.78))

# p_co2 <- ggplot() + geom_line(data=model_pco2_temp, aes(x=Year, y=`pCO2 (ppm)`, group=Group, color=Group), linewidth=1.2) + 
#   geom_point(data=hist_co2, aes(x=Year, y=CO2_ppm), size=2, color=pal[4], alpha=.5) +
#   scale_color_manual(values=c("#FEEDDE", "#FDBE85", "#FD8D3C", "#D94701", "grey")) +
#   labs(color="Global warming by 2100\nrelative to 1760 (°C)", x='', y="pCO2 (ppm)") +
#   theme(legend.position=c(0.3, 0.78))

# p <- p_temp +  p_co2 + plot_annotation(tag_levels = 'A')

fig_futureco2 <- ggplot() + 
  geom_line(data=model_pco2_temp %>% filter(Group!="PI"), aes(x=Year, y=`pCO2 (ppm)`, group=Group, color=Group)) + # model co2
  scale_color_manual(values=model_pal[-1])+
  geom_line(data=hist_co2, aes(x=Year, y=CO2_ppm), alpha=.5, color=pal[1])  + # observational co2
  labs(y="ppm", x="") + 
  annotate("rect", xmin = 1750, xmax=2020, ymin=280, ymax=410,
           fill="white", color="black", linetype="dotted", alpha=0)+
  theme(legend.position=c(0.4, 0.85), legend.direction = "horizontal", legend.background = element_blank(), legend.title = element_blank())

p_amoc <- ggplot() + geom_line(data=model_amoc, aes(x=Year, y=AMOC_max-model_amoc_baseline, group=Group, color=Group),
                               )+
  scale_color_manual(values=model_pal)+
  geom_ribbon(data=amoc_proxy, aes(x=Year,
                                   ymin=transport-amoc_proxy_baseline-SE,
                                   ymax=transport-amoc_proxy_baseline+SE,
                                   ), fill=pal[5], alpha=0.2) +
  geom_line(data=amoc_proxy, aes(x=Year, y=transport-amoc_proxy_baseline), color=pal[5]) +
  scale_x_date(limits = as.Date(c("1765-01-01","2100-01-01")))+
  labs(y="Sv", x="") + 
  theme(legend.position="none")

p_amoc_inset <- ggplot() + geom_line(data=amoc_proxy, aes(x=Year, y=transport-amoc_proxy_baseline), color=pal[5]) +
  geom_ribbon(data=amoc_proxy, aes(x=Year,
                                   ymin=transport-amoc_proxy_baseline-SE,
                                   ymax=transport-amoc_proxy_baseline+SE), fill=pal[5], alpha=0.2)+
  labs(x="", y="") + theme(panel.border = element_blank(),
                           panel.background = element_blank(),
                           axis.text = element_text(size=6))
p_amoc <- p_amoc + annotate("rect", xmin=as.Date("1765-01-01"), xmax=as.Date("2023-01-01"), ymin=-Inf, ymax=Inf,
                            fill="grey", alpha=.05)
p_amoc_inset <- p_amoc_inset+ annotate("rect", xmin = as.Date("1765-01-01"), xmax = as.Date("2000-01-01"), ymin=-7, ymax=0,
                       fill="white", color="black", alpha=.2) + theme(plot.margin=unit(c(0,0,0,0), "null"))

p_amoc <- p_amoc + inset_element(p_amoc_inset, left = 0.2, bottom = 0.2, right = 0.55, top = 0.55, align_to = 'full')


fig_dust_loading <- ggplot(data=gbl_dust_loading, aes(x=year)) + 
  geom_line(aes(y=global_dust_loading),color=pal[6]) + 
  geom_ribbon(aes(ymin=global_dust_loading-sigma, 
                  ymax=global_dust_loading+sigma), 
              alpha=0.2, fill=pal[6]) + 
  labs(y="Tg", x="Year")+
  scale_x_date(limits = as.Date(c("1765-01-01","2100-01-01")))

fig_dust_loading <- fig_dust_loading + annotate("rect", xmin=as.Date("1765-01-01"), xmax=as.Date("2023-01-01"), ymin=-Inf, ymax=Inf,
                            fill="grey", alpha=.05)

fig_futureco2 <- fig_futureco2+ annotate("rect",xmin=1765, xmax=2023, ymin=-Inf, ymax=Inf,
                                         fill="grey", alpha=.05)

p <- add_chronology(fig_co2) + fig_futureco2 +  add_chronology(fig_amoc) + p_amoc + 
  add_chronology(fig_iron) + fig_dust_loading+plot_layout(ncol = 2)

ggsave("output/divergence_deglacial_modern.pdf", width = 7, height = 9, units="in")
