library(tidyverse)
library(patchwork)
library(ggpubr)
library(rsvg) # convert to pdf

df_modern <- read_csv("data/modern_4pca.csv") %>% select(-1) %>% mutate(age  = "modern") %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
df_lgm <- read_csv("data/lgm_4pca.csv") %>% select(-1)  %>% mutate(age  = "LGM")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
df_future <- read_csv("data/future4_4pca.csv") %>% select(-1) %>% mutate(age  = "future")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
db <- rbind(df_modern, df_lgm, df_future) %>% pivot_longer(cols=bn:ss, names_to = "foram", values_to="biomass")
remove("df_modern", "df_lgm", "df_future")

db_longer <- db %>% mutate(biomass = log(biomass), fe=log(fe), po4=log(po4)) %>%
  pivot_longer(sst:stratification, names_to = "environment_variable",
                                 values_to = "environment_value")

# optimal niche
p_biomass_curve <- db_longer %>% filter(biomass > -4.5) %>% 
  ggplot(aes(x=environment_value, y=biomass, color=age,fill=age, group=age)) +
  facet_grid(foram~environment_variable, scales = "free_x")+ ylim(-4.5, -2)+
  geom_point(alpha=0.1)+
  geom_smooth(method="auto",se=F, fullrange=F)+
  ylab(expression(Log[10]~"biomass ("*Log[10]~"mmol C"~m^-3*")"))+ 
  scale_color_manual(values=c("#3D3B25FF", "#C71000FF", "#008EA0FF")) + 
  scale_fill_manual(values=c("#3D3B25FF", "#C71000FF", "#008EA0FF")) + 
  theme_bw() + 
  theme(text=element_text(family="Fira Sans"))+
  theme(strip.background = element_blank())

ggsave("output/env_biomass_curve.svg", p_biomass_curve, dpi=300, width=12, height=8)
rsvg::rsvg_pdf(svg="output/env_biomass_curve.svg",
               file ="output/env_biomass_curve.pdf")
file.remove("output/env_biomass_curve.svg")

## optimal niche, using density plot
p1 <- db_longer %>% filter(biomass > -4.5) %>% filter(foram=="bn") %>%
  ggdensity(x='environment_value', fill='age',facet.by = "environment_variable", nrow=1, scale="free", xlab="")+
  scale_fill_manual(values=c("#3D3B25FF", "#C71000FF", "#008EA0FF")) + theme(strip.background = element_blank())+
  theme(text=element_text(family="Fira Sans", size=10))

p2 <- db_longer %>% filter(biomass > -4.5) %>% filter(foram=="bs") %>%
  ggdensity(x='environment_value', fill='age',facet.by = "environment_variable", scale="free", nrow=1, legend="", xlab = "")+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  scale_fill_manual(values=c("#3D3B25FF", "#C71000FF", "#008EA0FF"))+
  theme(text=element_text(family="Fira Sans", size=10))

p3 <- db_longer %>% filter(biomass > -4.5) %>% filter(foram=="sn") %>%
  ggdensity(x='environment_value', fill='age',facet.by = "environment_variable", scale="free", nrow=1, legend="", xlab="")+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  scale_fill_manual(values=c("#3D3B25FF", "#C71000FF", "#008EA0FF"))+
  theme(text=element_text(family="Fira Sans", size=10))

p4 <- db_longer %>% filter(biomass > -4.5) %>% filter(foram=="ss") %>%
  ggdensity(x='environment_value', fill='age',facet.by = "environment_variable", scale="free", nrow=1, legend="")+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  scale_fill_manual(values=c("#3D3B25FF", "#C71000FF", "#008EA0FF"))+
  theme(text=element_text(family="Fira Sans", size=10))


p_density <- p1+p2+p3+p4 + plot_layout(nrow = 4)

ggsave("output/env_density_plot.svg", p_density, dpi=300, width=15, height=12)
rsvg::rsvg_pdf(svg="output/env_density_plot.svg",
               file ="output/env_density_plot.pdf")
file.remove("output/env_density_plot.svg")
