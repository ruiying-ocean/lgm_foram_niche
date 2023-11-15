## This script is to analyse the species optimal temperature change
## and its behind driver, here as the symbiont and spine trait

## Contact: rui.ying@bristol.ac.uk

## load visualisation packages
source("code/lib.R")

## read presaved Rdata
load("data/obs_smooth.Rdata")

Topt_diff <- obs_sp_smooth %>% thermal_opt(long_format=F) %>% mutate(diff = PI-LGM)

## export to csv
write_csv(Topt_diff, "data/Topt_sp_lgm.csv")

## merge with species trait data table
foram_sp_db <- read_csv("https://raw.githubusercontent.com/ruiying-ocean/lgm_foram_census/main/fg/foram_taxonomy.csv")

trait_info <-foram_sp_db %>%
   mutate(sp = map_vec(`Species name`, species_abbrev)) %>% select(sp, Symbiosis, Spine)

## check all species are included
common_sp <- intersect(Topt_diff$species, trait_info$sp)
if (length(common_sp) == length(Topt_diff$species)) {
  print("all species included")
}

Topt_diff <- Topt_diff %>% rename(sp=species) %>% left_join(trait_info, by = "sp") %>% distinct()

## anova analysis using symbiont and spine trait
Topt_diff <- Topt_diff %>% dplyr::filter(Symbiosis != 'undetermined',  Spine != 'underdetermined')


## merge symbiont obligate and symbiont bearing groups
## Topt_diff <- Topt_diff %>% mutate(Symbiosis = ifelse(Symbiosis == 'symbiont-obligate', 'symbiont-bearing', Symbiosis))

mod <- aov(diff ~ Spine + Symbiosis, data = Topt_diff)
summary(mod)

## comparison between groups
TukeyHSD(mod)

## ggpubr plot boxplot
library(ggpubr)

Topt_diff %>% mutate(ecogroup=paste(Symbiosis, Spine)) %>%
    ggboxplot(, x='ecogroup',y='diff',fill='ecogroup',palette = "jco",add='jitter')+
  stat_compare_means(method = "anova",label = "p.signif") +
  theme_publication()+
  theme(legend.position = "none", axis.title.x = element_blank())+
  ylab(expression(paste("PI - LGM optimal temperature (",degree,"C)")))

## save the figure
ggsave("output/figs2.png",width=8,height=4)
