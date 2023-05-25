## Alternative way to plot growth curve
library(patchwork)
library(ggpubr)
source("code/adaptation_data.R")
library(tidyverse)
library(quantregGrowth)

## fit model for each species
quantreg_model <- function(data, percentile_lvl=0.9){
  o <-gcrq(count ~ ps(sst), tau=percentile_lvl, data=data)
  return(o)
}

## get predicted values
predict_data <- function(data,tau=0.9,resolution=0.1){
  qr_model <- quantreg_model(data,percentile_lvl = tau)
  sp <- data$species %>% unique()
  age <- data$age %>% unique()
  
  min_sst <- min(data$sst)
  max_sst <- max(data$sst)
  sst <- seq(min_sst,max_sst,by=resolution)
  count <- charts(qr_model, k=sst, dataframe = T) %>% as.vector()
  chart <- tibble(count) %>% mutate(species = sp, age=age, sst=sst)
  return (chart)
}

## model each species in each age
## and save predicted values
data_list = list()
sp_list <- niche_data$species %>% unique()
age_list <- niche_data$age %>% unique()
i <- 1
for (sp_i in sp_list) {
  for (age_j in age_list) {
    tmp_raw <- niche_data %>% filter(age == age_j, species==sp_i) 
    tmp_predict <- predict_data(tmp_raw)
    data_list[[i]]<-tmp_predict
    i <- i + 1
    remove(tmp_raw)
    remove(tmp_predict)
  }
}

quantreg_pred_obs <- do.call("rbind", data_list)

## Do the same for GENIE output
data_list = list()
niche_genie <- niche_genie %>% rename(count=biomass) ## just for convenience
sp_list <- niche_genie$species %>% unique()
age_list <- niche_genie$age %>% unique()
i <- 1
for (sp_i in sp_list) {
  for (age_j in age_list) {
    tmp_raw <- niche_genie %>% filter(age == age_j, species==sp_i) 
    tmp_predict <- predict_data(tmp_raw, tau=0.9)
    data_list[[i]]<-tmp_predict
    i <- i + 1
    remove(tmp_raw)
    remove(tmp_predict)
  }
}

quantreg_pred_genie <- do.call("rbind", data_list)

## group species into functional types
# species order as functional groups
# symbiont-barren non-spinose
g1 <- c("N. pachyderma", "N. incompta", "T. quinqueloba")
# symbiont-barren spinose
g2 <- c("G. bulloides")
# symbiont-falcultative
g3 <- c("G. inflata", "N. dutertrei", "G. glutinata","G. menardii")
# symbiont-obligate
g4 <- c("G. ruber p", "G. ruber w", "T. sacculifer", "O. universa")

quantreg_pred_obs <- quantreg_pred_obs %>% mutate(function_group = case_when(
  species %in% g1 ~ "symbiont-barren non-spinose",
  species %in% g2 ~ "symbiont-barren spinose",
  species %in% g3 ~ "symbiont-falcultative non-spinose",
  species %in% g4 ~ "symbiont-obligate spinose",
), .after = age)


ggplot()  + geom_line(data=quantreg_pred_obs, aes(x=sst, y =count, color=age), linewidth=1)+
  facet_wrap(~species, scale="free_y") 

## plot
## raw data
p_s1 <- niche_data %>% ggplot(aes(x=sst, y=count))+  geom_point(aes(color=age), alpha=0.5, size=0.05)+
  facet_wrap(~species, scale="free_y")

## predicted data
p1 <- ggplot()  + geom_line(data=quantreg_pred_obs, aes(x=sst, y =count, color=age), linewidth=1)+
  facet_wrap(~species, scale="free_y") 

## p1 modify color and label
p1 <- p1 +  scale_color_manual(values=c("#0072B5FF", "#E18727FF", "#7876B1FF"),
                               labels=c("LGM","PI", "Modern")) + labs(x="Sea surface temperature (°C)", y="Abundance (#)")
p1 <- p1+ theme_bw() + theme(strip.background = element_blank(), strip.text = element_text(face = "italic"), legend.position = "none")

## statistic test

quantreg_pred_obs <- quantreg_pred_obs %>% mutate(coarse_sst=round(sst,0))

p2a <- quantreg_pred_obs %>% group_by(function_group, age, coarse_sst) %>% summarise(count = sum(count))%>% 
  ggplot() + geom_line(aes(x=coarse_sst, y=count, color=age)) + facet_wrap(~function_group)

p2b <- ggplot() +  geom_line(data=quantreg_pred_genie, aes(x=sst, y =count, color=age), linewidth=1)+
  facet_wrap(~species, scale="free_y")

p2a+p2b
