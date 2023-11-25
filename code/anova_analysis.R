## This script is to analyse the species optimal temperature change
## and its behind driver, here as the symbiont and spine trait

## Contact: rui.ying@bristol.ac.uk

## load visualisation packages
source("code/lib.R")

## read presaved Rdata
load("data/obs_smooth.Rdata")
load("data/genie_fg_smooth.Rdata")

Topt_sp <- obs_sp_r_smooth %>% thermal_opt(long_format=F,Topt_coef = 0.8) %>%
  mutate(Topt_mean_diff=PI_Topt_mean-LGM_Topt_mean)

Topt_fg <- obs_fg_r_smooth %>% thermal_opt(long_format=F,Topt_coef = 0.8) %>%
  mutate(Topt_mean_diff=PI_Topt_mean-LGM_Topt_mean)

modelled_Topt_fg <- genie_fg_smooth %>% filter(age=="lgm" | age=="pi") %>%
  thermal_opt(long_format=F,Topt_coef = 0.8)

## export to csv
write_csv(Topt_sp, "data/Topt_sp_lgm.csv")

## merge with species trait data table
foram_sp_db <- read_csv("https://raw.githubusercontent.com/ruiying-ocean/lgm_foram_census/main/fg/foram_taxonomy.csv")

trait_info <-foram_sp_db %>%
   mutate(sp = map_vec(`Species name`, species_abbrev)) %>% select(sp, Symbiosis, Spine)

## check all species are included
common_sp <- intersect(Topt_sp$species, trait_info$sp)
if (length(common_sp) == length(Topt_sp$species)) {
  print("all species included")
}

Topt_sp <- Topt_sp %>% rename(sp=species) %>% left_join(trait_info, by = "sp") %>% distinct()

## anova analysis using symbiont and spine trait
Topt_sp <- Topt_sp %>% dplyr::filter(Symbiosis != 'undetermined',  Spine != 'underdetermined') %>%
  dplyr::filter(Symbiosis != 'facultative',  Spine != 'underdetermined')

## merge symbiont obligate and symbiont bearing groups
mod <- aov(Topt_mean_diff ~ Spine + Symbiosis, data = Topt_sp)
summary(mod)

## ggpubr plot boxplot
library(ggpubr)

Topt_sp %>% mutate(ecogroup=paste(Symbiosis, Spine)) %>%
  filter(ecogroup!="symbiont-facultative spinose") %>%
  ggboxplot(x='ecogroup',y='Topt_mean_diff',fill='ecogroup',palette = "jco",add='jitter')+
  stat_compare_means(method = "anova",label = "p.signif") +
  theme_publication(base_size = 12)+
  theme(legend.position = "right", axis.text.x = element_blank())+
  ylab(expression(paste("∆ species optimal temperature (°C)"))) +
  xlab("Ecological group")

## save the figure to svg
ggsave(file = "output/figs5.svg", dpi = 300, width = 8, height = 5)
## convert to pdf
system("inkscape output/figs5.svg --export-pdf=output/figs5.pdf")
## rsvg solution
##system("rsvg-convert -f pdf output/fig1.svg > output/fig1.pdf")
## remove svg file
system("rm output/figs5.svg")
