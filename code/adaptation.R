library(tidyverse)
library(zoo)
select_species <- c("Orbulina universa", "Globigerina bulloides", "Neogloboquadrina pachyderma",
  "Neogloboquadrina dutertrei", "Neogloboquadrina incompta", "Globorotalia inflata", "Globigerinita glutinata",
  "Globigerinoides ruber w", "Globigerinoides ruber p", "Globorotalia menardii" ,
  "Turborotalita quinqueloba","Trilobatus sacculifer"
  )

lgm_climap <- read_csv("data/lgm_foram/raw/LGM_CLIMAP_wsst.csv")
lgm_atlantic <- read_csv("data/lgm_foram/raw/LGM_MARGO/LGM_MARGO_SouthAtlantic_count_wsst.csv")
lgm_pacific <- read_csv("data/lgm_foram/raw/LGM_MARGO/LGM_MARGO_pacific_count_wsst.csv")

names(lgm_climap) <- gsub(" [#]", "", names(lgm_climap), fixed=T)
names(lgm_climap) <- gsub(" [m]", "", names(lgm_climap), fixed=T)

#`G. truncatulinoides` is strangly high          
lgm_climap <- lgm_climap %>%
  select(SST,SSS,`O. universa`, `G. bulloides`,`G. inflata`,`G. quinqueloba`,
         `N. dutertrei`, `G. ruber w`, `G. ruber p`,
         `N. pachyderma d`, `N. pachyderma s`, `G. glutinata`,
         `G. menardii`,
         `G. sacculifer total`
         ) %>%
  rename(`N. incompta` = `N. pachyderma d`,
         `N. pachyderma` = `N. pachyderma s`,
         `T. quinqueloba`=`G. quinqueloba`,
         `T. sacculifer` = `G. sacculifer total`)
  
lgm_atlantic <- lgm_atlantic %>% rename(`Neogloboquadrina pachyderma` = `Neogloboquadrina pachyderma L`,
                        `Neogloboquadrina incompta`= `Neogloboquadrina pachyderma R`,
                        `Globigerinita glutinata` = `Globigerinita glutinata`,
                        `Globigerinoides ruber p` = `Globigerinoides ruber (pink)`,
                        `Globigerinoides ruber w` = `Globigerinoides ruber (white)`,
                        `Trilobatus sacculifer`  = `Globigerinoides sacc total`) %>%
  select(SST,SSS, all_of(select_species))

lgm_pacific <- lgm_pacific %>%
  rename(`Neogloboquadrina pachyderma` = `Neogloboquadrina pachyderma L`,
         `Neogloboquadrina incompta`= `Neogloboquadrina pachyderma R`,
         `Globigerinita glutinata` = `Globigerinita glutinata`,) %>%
  select(SST,SSS, any_of(select_species))

### ruber pink and white!
species_abbrev <- function(full_name, split_str = " ", sep_string=". "){
  splited_string <- str_split(full_name, split_str)[[1]]
  
  genus_name <- splited_string[1]
  sp_name <- splited_string[2]
  
  genus_abbrev <- str_sub(genus_name, 1, 1)
  combine_name <- paste(genus_abbrev, sp_name, sep = sep_string)
  
  if (length(splited_string) > 2){
    others <- str_split(full_name, split_str)[[1]][3]
    all_name <- paste0(combine_name, " ", others)
    return(all_name)
  } else{
    return(combine_name)
  }
}

lgm_pacific <- lgm_pacific %>% pivot_longer(cols =  !c(SST,SSS), names_to = "Species", values_to = "Count") %>%
  mutate(Species=sapply(Species, species_abbrev))

lgm_atlantic <- lgm_atlantic %>% pivot_longer(cols = !c(SST,SSS), names_to = "Species", values_to = "Count") %>%
  mutate(Species=sapply(Species, species_abbrev))

lgm_climap <- lgm_climap %>% pivot_longer(cols =  !c(SST,SSS), names_to = "Species", values_to = "Count") %>%
  mutate(Species=sapply(Species, species_abbrev))

lgm_niche <- rbind(lgm_pacific, lgm_atlantic, lgm_climap)
lgm_niche$age = "lgm"

### ForCenS
pi <- read_csv("data/modern_foram/raw/forcens_raw_count_wsst.csv")
pi_niche <-pi %>% select(c("Orbulina_universa", "Globigerina_bulloides", "Neogloboquadrina_pachyderma",
                    "Neogloboquadrina_dutertrei", "Neogloboquadrina_incompta", "Globoconella_inflata",
                    "Globigerinita_glutinata", "Turborotalita_quinqueloba",
                    "Trilobatus_sacculifer", "Globorotalia_menardii",
                    "Globigerinoides_white", "Globigerinoides_ruber"), SST, SSS) %>%
  rename(Globigerinoides_ruber_w=Globigerinoides_white,
         Globigerinoides_ruber_p=Globigerinoides_ruber)
pi_niche <- pi_niche %>%
  pivot_longer(cols = !c(SST,SSS), names_to = "Species", values_to = "Count") %>%
  mutate(Species=sapply(Species, species_abbrev, split_str="_"))

pi_niche$age = "pi"
#source("code/lombard_2009_model.R")
all_niche <- rbind(pi_niche,lgm_niche)

# n controls resolution
get_breaks <- function(v,n){
    min.v <- v %>% min(na.rm=T)
    max.v <- v %>% max(na.rm=T)
    step = (max.v - min.v)/(n-1)
    breaks = seq(min.v, max.v, step)
    return(breaks)
}

plot_niche <- function(n, quantile,variable){
  var = all_niche %>% pull(variable)
  b = get_breaks(var,n=n)#break point
  l <- rollmean(b,2) #midpoint of interval
  bin_df <- data.frame(l) %>% rowid_to_column("bin_id")
  tmp <- all_niche %>% mutate(bin_id = cut(all_niche[[variable]], breaks=b, labels=F)) %>%
    left_join(bin_df)
  
  # get highest values
  tmp <-  tmp %>% group_by(Species, age, bin_id) %>%
    summarise(Count, l, cutoff = quantile(Count, quantile, na.rm=T)) %>%
    filter(Count > cutoff)
  
  p <- tmp %>% ggplot(aes(x=l, y=Count, group=age)) + 
    stat_summary(fun=mean, geom="line", aes(color=age),alpha=0.25)+
    geom_smooth(aes(color=age),se=F)+
    facet_wrap(~Species, scales = "free_y")+
    scale_color_manual(values=c("#008EA0FF", "#C71000FF"),
                       labels=c("LGM","Modern"))
  return(p)
}

plot_niche(n=60, quantile=0.85, "SST") + xlab("Sea Surface Temperature (°C)") + 
  ylab("Abundance (#)")+
  theme_bw()+
  theme(text=element_text(family="Fira Sans"), strip.text = element_text(face="italic"))
ggsave("output/Niche_SST_adaptation.png", dpi=300)

plot_niche(n=120,  quantile=0.8,"SSS")+ xlab("Sea Suraface Salinity (PSU)") +
  ylab("Abundance (#)")+
  theme_bw()+
  theme(text=element_text(family="Fira Sans"), strip.text = element_text(face="italic"))
ggsave("output/Niche_SSS_adaptation.png", dpi=300)



