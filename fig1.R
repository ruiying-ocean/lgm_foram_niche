library(patchwork)
library(ggpubr)
library(tidyverse)
library(quantregGrowth)

lgm_fg <- read_csv("~/foram_core/tidy/lgm_fg_a_wsst.csv") %>% mutate(age="LGM")

lgm_fg <- lgm_fg %>% pivot_longer(cols = `Symbiont-barren Non-Spinose`:`Symbiont-obligate Spinose`,
                        names_to = "Group",values_to = "Abundance")

pi_fg <- read_csv("~/foram_core/tidy/forcens_fg_a_wsst.csv") %>% mutate(age="PI")
pi_fg <- pi_fg %>% pivot_longer(cols = `Symbiont-obligate Spinose`:`Symbiont-barren Spinose`,
                                  names_to = "Group",values_to = "Abundance")

all_fg <- rbind(pi_fg, lgm_fg)

## fit model for each species
quantreg_model <- function(data, percentile_lvl=0.9){
  o <-gcrq(Abundance ~ ps(SST), tau=percentile_lvl, data=data)
  return(o)
}

## get predicted values
predict_data <- function(data,tau=0.95,resolution=0.1){
  qr_model <- quantreg_model(data,percentile_lvl = tau)
  sp <- data$Group %>% unique()
  age <- data$age %>% unique()
  
  min_sst <- min(data$SST, na.rm=T)
  max_sst <- max(data$SST, na.rm = T)

  sst <- seq(min_sst,max_sst,by=resolution)
  count <- charts(qr_model, k=sst, dataframe = T) %>% as.vector()
  chart <- tibble(count) %>% mutate(species = sp, age=age, sst=sst)
  return (chart)
}

## model each species in each age
## and save predicted values
data_list = list()
gr_list <- all_fg$Group %>% unique()
age_list <- all_fg$age %>% unique()
i <- 1
for (gr_i in gr_list) {
  for (age_j in age_list) {
    tmp_raw <- all_fg %>% filter(age == age_j, Group==gr_i) 
    tmp_predict <- predict_data(tmp_raw,resolution = 0.5)
    data_list[[i]]<-tmp_predict
    i <- i + 1
    remove(tmp_raw)
    remove(tmp_predict)
  }
}


quantreg_pred_obs <- do.call("rbind", data_list)
all_fg <- all_fg %>% rename(species=Group)
# quantreg_pred_obs <- quantreg_pred_obs %>% filter(species!="Symbiont-facultative Non-Spinose")
# all_fg<-all_fg %>% filter(species!="Symbiont-facultative Non-Spinose") 

p1a <- ggplot()+
  geom_point(aes(x=SST, y=Abundance, color=age), data=all_fg, size=0.1, alpha=0.5) +
  geom_line(data=quantreg_pred_obs,aes(x=sst, y=count, color=age),linewidth=1) + 
  facet_wrap(~species,scales = "free_y",nrow=1)

p1a <- p1a +  scale_color_manual(values=c("#0072B5FF", "#DC143C", "#7876B1FF"),
                            labels=c("LGM","Holocene", "Modern")) + labs(x="", y="Abundance (#)")
p1a <- p1a+ theme_bw() + theme(strip.background = element_blank(), strip.text = element_text(face = "italic"), legend.position = "none")
p1a
