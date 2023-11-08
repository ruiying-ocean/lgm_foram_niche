## This script is to analyse the species optimal temperature change
## and its behind driver, here as the symbiont and spine trait

## Contact: rui.ying@bristol.ac.uk

## load visualisation packages
source("code/lib.R")

## read presaved Rdata
load("data/obs_smooth.Rdata")

Topt_diff <- obs_sp_smooth %>% thermal_opt(long_format=F) %>% mutate(diff = PI-LGM)

## merge with species trait data table
foram_sp_db <- read_csv("https://raw.githubusercontent.com/ruiying-ocean/lgm_foram_census/main/fg/foram_sp_db.csv")

trait_info <-foram_sp_db %>%
   mutate(sp = map_vec(Name, species_abbrev)) %>% select(sp, Symbiosis, Spinose)

## check all species are included
common_sp <- intersect(Topt_diff$species, trait_info$sp)
if (length(common_sp) == length(Topt_diff$species)) {
  print("all species included")
}

Topt_diff <- Topt_diff %>% rename(sp=species) %>% left_join(trait_info, by = "sp") %>% distinct()

Topt_diff <- Topt_diff %>% mutate(ecogroup = case_when(
  Symbiosis == "Yes" & Spinose == "Yes"~ "Symbiont-obligate Spinose",
  Symbiosis == "No" & Spinose == "Yes" ~ "Symbiont-barren Spinose",
  Symbiosis == "No" & Spinose == "No" ~ "Symbiont-barren Non-spinose",
  TRUE ~ "Symbiont-facultative Non-spinose")) %>%
  filter(ecogroup != "Symbiont-facultative Non-spinose")

## anova analysis using symbiont and spine trait
mod <- aov(diff ~ ecogroup, data = Topt_diff)
summary(mod)

## comparison between groups
TukeyHSD(mod)

## ggpubr plot boxplot
library(ggpubr)

ggboxplot(Topt_diff, x='ecogroup',y='diff',fill='ecogroup',palette = "jco",add='jitter')+
  stat_compare_means(method = "anova",label = "p.signif") +
  theme_publication()+
  theme(legend.position = "none", axis.title.x = element_blank())+
  ylab(expression(paste("PI - LGM optimal temperature (",degree,"C)")))

## save the figure
ggsave("output/figs10.png",width=8,height=4)
