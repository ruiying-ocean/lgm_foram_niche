library(tidyverse)
library(ggpubr)
library(ncdf4)

##  ----------- observational CO2/TEMP data -------------------
## historical temperature and CO2 data source: https://ourworldindata.org/co2-and-other-greenhouse-gas-emissions
hist_temp <- read_csv("model/model_drived/observed_global_temperature_anamoly.csv")
hist_co2 <- read_csv("model/model_drived/observed_co2.csv")
names(hist_temp)[2] <- "Global_Air_Temp_Anomaly"

## convert data baseline from 1961_1990 to 1850-1900
hist_mean_1960_1990 <- hist_temp %>%
  filter(Year > 1960 & Year < 1990) %>%
  summarise(mean_1960_1990 = mean(Global_Air_Temp_Anomaly)) %>%
  pull(mean_1960_1990)

hist_mean_1850_1900 <- hist_temp %>%
  filter(Year > 1850 & Year < 1900) %>%
  summarise(mean_1850_1900 = mean(Global_Air_Temp_Anomaly)) %>%
  pull(mean_1850_1900)

delta <- hist_mean_1960_1990 - hist_mean_1850_1900
hist_temp <- hist_temp %>% mutate(Global_Air_Temp_Anomaly = Global_Air_Temp_Anomaly + delta)

## --------------- cGENIE Model CO2/TEMP data -------------------
model_pco2_temp <- read_csv("model/model_drived/model_pCO2_temperature.csv")
model_pco2_temp$Year <- model_pco2_temp$Year - 0.5
model_pco2_temp <- model_pco2_temp %>% mutate(Group = recode(Group, "1p5" = "+1.5°C", "2" = "+2°C", "3" = "+3°C", "4" = "+4°C"))
model_pco2_temp$Group <- factor(model_pco2_temp$Group, levels = c("historical", "+1.5°C", "+2°C", "+3°C", "+4°C"))

## as the data, set average 1960-1990 as standard
genie_mean_1850_1900 <- model_pco2_temp %>%
  filter(Year > 1850 & Year < 1900, Group == "historical") %>%
  summarise(mean_1850_1900 = mean(Temperature)) %>%
  pull(mean_1850_1900)

model_pco2_temp <- model_pco2_temp %>%
  group_by(Group) %>%
  mutate(delta_T = Temperature - genie_mean_1850_1900) %>%
  ungroup()

## ---------- CMIP6 CO2 data -------------------
## cmip_co2 <- function(dir, origin_day = "0000-01-01") {
##   nc <- nc_open(dir)
##   # Read specific variables, dimensions, or attributes from the NetCDF file
##   co2 <- ncvar_get(nc, "mole_fraction_of_carbon_dioxide_in_air")

##   ## extract the global value
##   global_co2 <- co2[1, ]

##   # time since 0-1-1
##   time <- ncvar_get(nc, "time")

##   # Convert time to real date
##   time_origin <- as.POSIXct(origin_day, tz = "UTC")
##   ymd <- time_origin + as.difftime(time, unit = "days")

##   nc_close(nc)

##   return(data.frame(date = ymd, carbon_dioxide_concentration = global_co2))
## }

## ssp_co2_hist <- cmip_co2("~/Science/cmip_ts/mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gr1-GMNHSH_0000-2014.nc")
## ssp_co2_ssp126 <- cmip_co2("~/Science/cmip_ts/mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-IMAGE-ssp126-1-2-1_gr1-GMNHSH_2015-2500.nc", "1850-01-01 00:00:00")
## ssp_co2_ssp245 <- cmip_co2("~/Science/cmip_ts/mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-MESSAGE-GLOBIOM-ssp245-1-2-1_gr1-GMNHSH_2015-2500.nc", "1850-01-01 00:00:00")
## ssp_co2_ssp370 <- cmip_co2("~/Science/cmip_ts/mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-AIM-ssp370-1-2-1_gr1-GMNHSH_2015-2500.nc", "1850-01-01 00:00:00")
## ssp_co2_ssp585 <- cmip_co2("~/Science/cmip_ts/mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-REMIND-MAGPIE-ssp585-1-2-1_gr1-GMNHSH_2015-2500.nc", "1850-01-01 00:00:00")

## save to R data for reproducbility
## save(ssp_co2_hist, ssp_co2_ssp126, ssp_co2_ssp245, ssp_co2_ssp370, ssp_co2_ssp585, file = "data/ssp_co2.RData")

load("data/ssp_co2.RData")

ssp_co2_hist <- ssp_co2_hist %>% mutate(label = "historical")
ssp_co2_ssp126 <- ssp_co2_ssp126 %>% mutate(label = "ssp126")
ssp_co2_ssp245 <- ssp_co2_ssp245 %>% mutate(label = "ssp245")
ssp_co2_ssp370 <- ssp_co2_ssp370 %>% mutate(label = "ssp370")
ssp_co2_ssp585 <- ssp_co2_ssp585 %>% mutate(label = "ssp585")

ssp_co2_data <- rbind(ssp_co2_hist, ssp_co2_ssp126, ssp_co2_ssp245, ssp_co2_ssp370, ssp_co2_ssp585)
ssp_co2_data <- ssp_co2_data %>% filter(as.POSIXlt(date)$year + 1900 <= 2100, date >= as.Date("1850-01-01"))
ssp_co2_data <- ssp_co2_data %>% as_tibble()
## ---------- CMIP6 TEMP data -------------------

## ts_hist <- read_csv("~/Science/cmip_ts/ts_historical.csv")
## ts_ssp126 <- read_csv("~/Science/cmip_ts/ts_ssp126.csv")
## ts_ssp245 <- read_csv("~/Science/cmip_ts/ts_ssp245.csv")
## ts_ssp370 <- read_csv("~/Science/cmip_ts/ts_ssp370.csv")
## ts_ssp585 <- read_csv("~/Science/cmip_ts/ts_ssp585.csv")

## save as Rdata for reproducibility
## save(ts_hist, ts_ssp126, ts_ssp245, ts_ssp370, ts_ssp585, file = "data/ssp_temp.RData")

load("data/ssp_temp.RData")

## Combine all dataframes
cmip_ts <- rbind(ts_hist, ts_ssp126, ts_ssp245, ts_ssp370, ts_ssp585)

cmip_ts <- cmip_ts %>% rename(Monthly_Surf_Temp = Mean)

# Modify the Date column to remove the time of day
cmip_ts$Date <- as.Date(cmip_ts$Date)

# Extract month and year from Date column
cmip_ts$Month <- format(cmip_ts$Date, "%m")
cmip_ts$Year <- format(cmip_ts$Date, "%Y")

## calculate model average
cmip_ts_avg <- cmip_ts %>%
  group_by(Year, Month, experiment) %>%
  summarise(
    Mean = mean(Monthly_Surf_Temp, na.rm = TRUE),
    SD = sd(Monthly_Surf_Temp, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(Date = as.Date(paste(Year, Month, "15", sep = "-")))

# Calculate pre-industrial average
pre_industrial_avg <- cmip_ts_avg %>%
  filter(Year > 1850 & Year < 1900) %>%
  summarise(PreIndustrialAvg = mean(Mean, na.rm = TRUE)) %>%
  pull(PreIndustrialAvg)

cmip_ts_avg <- cmip_ts_avg %>% mutate(anomaly = Mean - pre_industrial_avg)


## --------------- Ploting starts here -------------------
## some ideas from https://rstudio-pubs-static.s3.amazonaws.com/284329_c7e660636fec4a42a09eed968dc47f32.html
theme_set(
  theme(
    text = element_text(family = "helvetica", size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray50", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(),
    legend.position = "bottom",
  )
)

## GENIE pCO2
p1 <- ggplot() +
  geom_line(data = model_pco2_temp, aes(x = Year, y = `pCO2 (ppm)`, group = Group, color = Group), linewidth = 1) +
  geom_point(data = hist_co2, aes(x = Year, y = CO2_ppm), size = 2, alpha = 0.5) +
  labs(
    color = "", x = "", y = "",
    title = expression("Global " ~ CO[2] ~ " concentration (ppm)")
  )

## GENIE CO2 concentration
p2 <- ssp_co2_data %>% ggplot(aes(x = date, y = carbon_dioxide_concentration)) +
  geom_line(aes(color = label), linewidth = 1)

p2 <- p2 + labs(
  color = "", x = "", y = "",
  title = expression("Global " ~ CO[2] ~ " concentration (ppm)")
)

## GENIE temperature relative to 1960-1990 average

p3 <- ggplot() +
  geom_line(data = model_pco2_temp, aes(x = Year, y = delta_T, group = Group, color = Group), linewidth = 1) +
  geom_line(data = hist_temp, aes(x = Year, y = Global_Air_Temp_Anomaly), linewidth = 1, alpha = 0.5) +
  labs(color = "", x = "", y = "", title = "Global surafce temperature anomaly (°C)")

## CMIP6 temperature anomaly

p4 <- ggplot(cmip_ts_avg, aes(x = Date, y = anomaly, group = experiment)) +
  geom_line(aes(color = experiment)) +
  geom_ribbon(aes(ymin = anomaly - SD, ymax = anomaly + SD, fill = experiment), alpha = 0.2) +
  labs(
    x = "", y = "", color = "", fill = "", title =
      "Global surface temperature anomaly (°C)"
  ) +
  guides(color = "none")

# p3 <- p3

## make x and y axis consistent
p3 <- p3 + xlim(1850, 2100) + scale_y_continuous(breaks = c(0, 1.5, 2, 3, 4), limits = c(-2, 8))
p4 <- p4 + ylim(-2, 8)
p1 <- p1 + xlim(1850, 2100) + ylim(250, 1200)
p2 <- p2 + ylim(250, 1200)

p3 <- p3 + geom_text(aes(x = 2060, y = -0.5, label = "1850-1900 Baseline"), color = "grey", size = 4)
p4 <- p4 + geom_text(aes(x = as.Date("2060-01-01"), y = -0.5, label = "1850-1900 Baseline"), color = "grey", size = 4)

## label observed data
p1 <- p1 + geom_text(aes(x = 1910, y = 380, label = "CO2 Observation"), color = "grey", size = 4)
p3 <- p3 + geom_text(aes(x = 1920, y = 1, label = "Temperature Observation"), color = "grey", size = 4)

## label model at the upper right corner
p1 <- p1 + geom_text(aes(x = 2060, y = 900, label = "cGENIE"), color = "black", size = 4)
p2 <- p2 + geom_text(aes(x = as.POSIXct("2060-01-01"), y = 1050, label = "CMIP6"), color = "black", size = 4)
p3 <- p3 + geom_text(aes(x = 2060, y = 4, label = "cGENIE"), color = "black", size = 4)
p4 <- p4 + geom_text(aes(x = as.Date("2060-01-01"), y = 6, label = "CMIP6"), color = "black", size = 4)

## remove the legends of first two plots
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")

## arrange the plots
p <- ggpubr::ggarrange(p1, p2, p3, p4,
  ncol = 2, nrow = 2,
  labels = c("(a)", "(b)", "(c)", "(d)"), heights = c(4, 5),
  label.x = 0.12, label.y = 0.85
)

## saving figures
ggsave("output/future_temperature_co2.png", p, width = 10, height = 7, dpi = 300)
