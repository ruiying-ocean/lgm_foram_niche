library(tidyverse)
library(santoku)
source("code/adaptation_data.R")

bin_along_var <- function(data, var, n, ...){
  env_var = data %>% pull(var)
  data <- data %>% mutate(bin_lvl=chop_evenly(env_var, n))
  step <- (max(env_var) - min(env_var))/(n)
  bin_mid_pnt <- min(env_var) + step/2 * seq(1, 2*n, 2)
  bin_lvl <- tab_evenly(env_var, n) %>% names()
  bin_df <- data.frame(bin_mid_pnt, bin_lvl)
  data <- data %>% left_join(bin_df)
  
  ## take average for abundance in each bin
  data <- data %>% group_by(Species, age, bin_mid_pnt) %>%
    summarise(mean_abundance = mean(biomass, na.rm = T))
  
  return (data)
}

foram_long_names <- c(
  `bn` = "symbiont-barren\nnon-spinose",
  `bs` = "symbiont-barren\nspinose",
  `sn` = "symbiont-facultative\nnon-spinose",
  `ss` = "symbiont-obligate\nspinose"
)

font_size=14

p1 <- niche_genie %>% bin_along_var(., "sst", n=45) %>% ggplot(aes(x=bin_mid_pnt, y=mean_abundance)) + geom_point(aes(color=age),size=1)+
  geom_line(aes(color=age),linewidth=.85, alpha=0.8)+facet_wrap(~Species, scale="free_y", labeller = as_labeller(foram_long_names), nrow = 1) 

p1 <- p1 +  scale_color_manual(values=c("#0072B5FF", "#E18727FF", "#7876B1FF"),
                      labels=c("LGM", "PI")) + labs(x="Sea surface temperature (°C)", y="Biomass (mmol C/m3)")

p1 <- p1 + theme_bw(base_size = font_size) + theme(strip.background = element_blank(), strip.text = element_text(face = "italic"))

p2 <- niche_genie %>% bin_along_var(., "chl_total", n=45) %>% ggplot(aes(x=bin_mid_pnt, y=mean_abundance)) + geom_point(aes(color=age),size=1)+
  geom_line(aes(color=age),linewidth=.85, alpha=0.8)+facet_wrap(~Species, scale="free_y", labeller = as_labeller(foram_long_names), nrow = 1) 

p2 <- p2 +  scale_color_manual(values=c("#0072B5FF", "#E18727FF", "#7876B1FF"),
                                 labels=c("LGM", "PI")) + labs(x="Total Chl (mg/m3)", y="Biomass (mmol C/m3)")
p2 <- p2 + theme_bw(base_size = font_size) + theme(strip.background = element_blank(), strip.text = element_blank())

p <-p1 / p2 & ylab(NULL)
ylab <- p1$labels$y
# Use the tag label as a y-axis label
p <- wrap_elements(p) +
  labs(tag = ylab) +
  theme(
    plot.tag = element_text(size = font_size, angle = 90),
    plot.tag.position = "left"
  )

p %>% ggsave("output/model_niche.jpg",., dpi=400, width=10, height=5)
niche_genie %>% ggplot(aes(x=sst,y=biomass)) + geom_point() + facet_wrap(~Species)
