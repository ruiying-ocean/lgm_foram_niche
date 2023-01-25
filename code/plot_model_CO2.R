library(tidyverse)
library(patchwork)
library(hrbrthemes)
library(RColorBrewer)

# historical graph
# temperature and CO2 data source: https://ourworldindata.org/co2-and-other-greenhouse-gas-emissions
past_temp <- read_csv("data/model_data/observed_global_temperature_anamoly.csv")
past_co2 <- read_csv("data/model_data/observed_co2.csv")
names(past_temp)[2] <- "Global_Air_Temp_Anomaly"

# model data
df <- read_csv("data/model_data/model_pCO2_temperature.csv")
#df <- df %>% mutate(Age_group = as.factor(if_else(Year < 2022, "History", "Future")))
df$Year <- df$Year - 0.5
df$Group <- as.factor(df$Group)

# as the data, set average 1960-1990 as standard
Average_standard <- df %>% filter(Year > 1960 & Year < 1990, Group=='PI') %>% summarise(avg_temp = mean(Temperature))
df <- df %>% group_by(Group) %>% mutate(delta_T= Temperature - Average_standard$avg_temp) %>% ungroup()

p_temp <- ggplot() + geom_line(data=df, aes(x=Year, y=delta_T, group=Group, color=Group), linewidth=1.2) + 
  geom_line(data=past_temp, aes(x=Year, y=Global_Air_Temp_Anomaly)) +
  theme_ipsum() + scale_color_ipsum() + 
  labs(color="Global warming by 2100\nrelative to 1760 (°C)", x='', y="Tempearature Relative to 1960-1990 (°C)") + 
  theme(legend.position=c(0.3, 0.78))

p_co2 <- ggplot() + geom_line(data=df, aes(x=Year, y=`pCO2 (ppm)`, group=Group, color=Group), linewidth=1.2) + 
  geom_point(data=past_co2, aes(x=Year, y=CO2_ppm), size=1.5) +
  theme_ipsum() + scale_color_ipsum()  + labs(color="Global warming by 2100\nrelative to 1760 (°C)", x='', y="pCO2 (ppm)") +
  theme(legend.position=c(0.3, 0.78))

p <- p_temp +  p_co2 + plot_annotation(tag_levels = 'A')

png("./output/temperature_gradient.jpg", res=300, width = 10, height = 5, unit="in")
print(p)
dev.off()
