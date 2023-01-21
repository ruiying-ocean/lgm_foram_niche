library(tidyverse)
library(factoextra)
library(patchwork)
library(geometry)

df_modern <- read_csv("data/modern_4pca.csv") %>% select(-1) %>% mutate(age  = "preindustrial") %>% mutate(age = factor(age, levels=c("future", "preindustrial", "LGM")))
df_lgm <- read_csv("data/lgm_4pca.csv") %>% select(-1)  %>% mutate(age  = "LGM")  %>% mutate(age = factor(age, levels=c("future", "preindustrial", "LGM")))
df_future <- read_csv("data/future4_4pca.csv") %>% select(-1) %>% mutate(age  = "future")  %>% mutate(age = factor(age, levels=c("future", "preindustrial", "LGM")))
db <- rbind(df_modern, df_lgm, df_future) %>% pivot_longer(cols=bn:ss, names_to = "foram", values_to="biomass")
remove("df_modern", "df_lgm", "df_future")

db <- db %>% mutate(biomass=log(biomass)) %>% filter(biomass>-6)

db %>% filter(foram=="bn") %>% ggplot(aes(x=sst, y=biomass, color=age, group=age)) +
  geom_point(alpha=0.2)+
  geom_smooth(method="gam", se=T, fullrange=F)

# niche changed!