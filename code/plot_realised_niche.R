library(tidyverse)
library(factoextra)
library(patchwork)
library(geometry)
## library(wesanderson)
library(pals)
library(ggsci)
library(egg)

df_modern <- read_csv("data/modern_4pca.csv") %>% select(-1) %>% mutate(age  = "modern") %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
df_lgm <- read_csv("data/lgm_4pca.csv") %>% select(-1)  %>% mutate(age  = "LGM")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
df_future <- read_csv("data/future4_4pca.csv") %>% select(-1) %>% mutate(age  = "future")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
db <- rbind(df_modern, df_lgm, df_future) %>% pivot_longer(cols=bn:ss, names_to = "foram", values_to="biomass")
remove("df_modern", "df_lgm", "df_future")

db <- db %>% mutate(biomass=log(biomass)) %>% filter(biomass>-9)

foram_names <- c(
  `bn` = "symbiont-barren non-spinose",
  `bs` = "symbiont-barren spinose",
  `sn` = "symbiont-facultative non-spinose",
  `ss` = "symbiont-obligate spinose"
)

p1 <- db %>% ggplot(aes(x=sst, y=biomass, color=age,fill=age, group=age)) +
  geom_point(alpha=0.3)+
  geom_smooth(method="loess",se=T, fullrange=F) +
  facet_wrap(~foram,  strip.position = "bottom", labeller = as_labeller(foram_names)) +
  xlab("SST (°C)")+ 
  ylab(expression(Log[10]~"biomass ("*Log[10]~"mmol C"~m^-3*")"))+ 
  scale_color_manual(values=c("#3D3B25FF", "#C71000FF", "#008EA0FF")) + 
  scale_fill_manual(values=c("#3D3B25FF", "#C71000FF", "#008EA0FF")) + 
  theme_bw() + 
  theme(text=element_text(family="Fira Sans"))+
  theme(strip.placement = "outside", strip.background = element_blank())

p2 <- db %>% ggplot(aes(x=prey, y=biomass, color=age,fill=age, group=age)) +
  geom_point(alpha=0.3)+
  geom_smooth(method="loess",se=T, fullrange=F) +
  facet_wrap(~foram,  strip.position = "bottom", labeller = as_labeller(foram_names)) +
  xlab(expression("Preferred prey biomass (mmol C"~m^-3*")"))+ 
  ylab(expression(Log[10]~"biomass ("*Log[10]~"mmol C"~m^-3*")"))+ 
  scale_color_manual(values=c("#3D3B25FF", "#C71000FF", "#008EA0FF")) + 
  scale_fill_manual(values=c("#3D3B25FF", "#C71000FF", "#008EA0FF")) + 
  theme_bw() + 
  theme(text=element_text(family="Fira Sans"))+
  theme(strip.placement = "outside", strip.background = element_blank())

p1 <- tag_facet(p1)
p2 <- tag_facet(p2)

ggsave("output/realised_thermal_niche.jpg", p1, dpi=300)
ggsave("output/realised_food_niche.jpg", p2, dpi=300)

