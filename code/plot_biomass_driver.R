library(tidyverse)
library(factoextra)
library(patchwork)
library(geometry)
## library(wesanderson)
library(pals)
library(ggsci)

df_modern <- read_csv("data/modern_4pca.csv") %>% select(-1) %>% mutate(age  = "modern") %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
df_lgm <- read_csv("data/lgm_4pca.csv") %>% select(-1)  %>% mutate(age  = "LGM")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
df_future <- read_csv("data/future4_4pca.csv") %>% select(-1) %>% mutate(age  = "future")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
db <- rbind(df_modern, df_lgm, df_future) %>% pivot_longer(cols=bn:ss, names_to = "foram", values_to="biomass")
remove("df_modern", "df_lgm", "df_future")

db <- db %>% mutate(biomass=log(biomass)) %>% filter(biomass>-6)

p1 <- db %>% ggplot(aes(x=sst, y=biomass, color=age,fill=age, group=age)) +
  geom_point(alpha=0.3)+
  geom_smooth(method="gam", se=T, fullrange=F) +
  facet_wrap(~foram) + theme_bw() + 
  theme(text=element_text(family="Fira Sans"))+
  scale_color_jama() +  scale_fill_jama()

p2 <- db %>% ggplot(aes(x=prey, y=biomass, color=age,fill=age, group=age)) +
  geom_point(alpha=0.3)+
  geom_smooth(method="gam", se=T, fullrange=F) +
  facet_wrap(~foram) + theme_bw() + 
  theme(text=element_text(family="Fira Sans"))+
  scale_color_jama() +  scale_fill_jama()

ggsave("output/thermal_niche.jpg", p1, dpi=300)
ggsave("output/food_niche.jpg", p2, dpi=300)

db %>% ggplot(aes(x=log(po4), y=biomass, color=age,fill=age, group=age)) +
  geom_point(alpha=0.3)+
  geom_smooth(method="gam", se=T, fullrange=F) +
  facet_wrap(~foram) + theme_bw() + 
  theme(text=element_text(family="Fira Sans"))+
  scale_color_jama() +  scale_fill_jama()

db %>% ggplot(aes(x=lat, y=biomass, color=age,fill=age, group=age)) +
  geom_point(alpha=0.3)+
  geom_smooth(se=T, fullrange=F) +
  facet_wrap(~foram) + theme_bw() + 
  theme(text=element_text(family="Fira Sans"))+
  scale_color_jama() +  scale_fill_jama()
