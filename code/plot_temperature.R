library(tidyverse)
library(patchwork)
library(hrbrthemes)
library(RColorBrewer)

df <- read_csv("data/Future_temp.csv")
df$Group <- as.factor(df$Group)
df$pCO2 <- df$pCO2 * 1E6
p2 <- ggplot(data=df, aes(x=Year)) + geom_line(aes(y=Temperature, group=Group, color=Group), linewidth=1.2) + ylab("Tempearature (C)") +theme_ipsum() + scale_color_ipsum() + theme(legend.position="none")
p1 <- ggplot(data=df, aes(x=Year)) + geom_line(aes(y=pCO2, group=Group, color=Group), linewidth=1.2) + ylab("pCO2 (ppm)")  +theme_ipsum() + scale_color_ipsum()
p <- p1 + p2 + plot_annotation(tag_levels = 'A')


# historical graph
past_temp <- read_csv("data/historical_global_temperature_anamoly.csv")
past_co2 <- read_csv("data/historical_co2.csv")

my_colors <- brewer.pal(5, "Reds")
my_colors <- my_colors[-1]

zip_colors <- c("1.5" = my_colors[1],
            "2" = my_colors[2],
            "3" = my_colors[3],
            "4" = my_colors[4])

p_co2 <- ggplot(data=past_co2, aes(x=Year)) +  geom_point(aes(y=CO2_ppm), alpha=0.2)
p_co2 + geom_segment(aes(x = 1765, y = 278, xend = 1960, yend = 316)) + 
  geom_segment(aes(x= 1960, y = 316, xend=2022, yend=420)) +
  geom_segment(aes(x= 2022, y = 420, xend=2100, yend=420, color="1.5")) +
  geom_segment(aes(x= 2022, y = 420, xend=2100, yend=490, color="2")) +
  geom_segment(aes(x= 2022, y = 420, xend=2100, yend=650, color="3")) +
  geom_segment(aes(x= 2022, y = 420, xend=2100, yend=833, color="4")) +
  geom_text(aes(x=Inf, y=420, label="420"), hjust=1.5, vjust=-1, size=2.5) +
  geom_text(aes(x=Inf, y=490, label="490"), hjust=1.5, vjust=-1, size=2.5) +
  geom_text(aes(x=Inf, y=650, label="650"), hjust=1.5, vjust=-1, size=2.5) +
  geom_text(aes(x=Inf, y=833, label="833"), hjust=1.5, vjust=-1, size=2.5) +
  scale_color_manual(values = zip_colors) + 
  labs(x = "Year",
       y = "CO2 (ppm)",
       color = "") +  theme_ipsum() + theme(legend.position = c(0.1, 0.9))

p_temp <- ggplot(data=past_temp, aes(x=Year)) +  geom_line(aes(y=Temperature_Anomaly))

#png("./output/temperature_gradient.jpg", res=300, width = 8, height = 4, unit="in")
#print(p)
#dev.off()
