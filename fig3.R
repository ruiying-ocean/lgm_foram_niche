library(patchwork)
library(ggpubr)
library(tidyverse)
library(quantregGrowth)

lgm_sp <- read_csv("~/foram_core/tidy/lgm_sp_a_wsst.csv") %>% mutate(age="LGM")

lgm_sp <- lgm_sp %>% pivot_longer(cols = `O. universa`:`G. theyeri`,
                                  names_to = "species",values_to = "Abundance")

pi_sp <- read_csv("~/foram_core/tidy/forcens_sp_a_wsst.csv") %>% mutate(age="PI")


pi_sp <- pi_sp %>% pivot_longer(cols = `D. anfracta`:`H. digitata`,
                                names_to = "species",values_to = "Abundance")

top_lgm_sp <- lgm_sp %>% group_by(species) %>% summarise(ab = mean(Abundance, na.rm=T))%>%
  arrange(desc(ab)) %>% head(30) %>% pull(species)

top_pi_sp <- pi_sp %>% group_by(species) %>% summarise(ab = mean(Abundance, na.rm=T))%>%
  arrange(desc(ab)) %>% head(30) %>% pull(species)

# select_species <- intersect(top_lgm_sp, top_pi_sp) 


select_species <- c("G. bulloides", "N. pachyderma",
                    "N. dutertrei", "N. incompta", "G. inflata", "G. glutinata",
                    "G. ruber albus", "G. ruber ruber","G. truncatulinoides",
                    "P. obliquiloculata", "T. quinqueloba","T. sacculifer")


subset_columns <- c("SST", "species", "Abundance", "age")
all_sp <- rbind(pi_sp[subset_columns], lgm_sp[subset_columns])
all_sp <- all_sp %>% filter(species %in% select_species)

## fit model for each species
quantreg_model <- function(data, percentile_lvl=0.9){
  o <-gcrq(Abundance ~ ps(SST), tau=percentile_lvl, data=data)
  return(o)
}

## get predicted values
predict_data <- function(data,tau=0.95,resolution=0.1){
  
  data <- data %>% drop_na(SST)
  
  qr_model <- quantreg_model(data,percentile_lvl = tau)
  
  sp <- data$species %>% unique()
  age <- data$age %>% unique()

  print(sp)
  print(age)
  min_sst <- min(data$SST, na.rm=T)
  max_sst <- max(data$SST, na.rm = T)
  print(min_sst)
  print(max_sst)

  sst <- seq(min_sst,max_sst,by=resolution)
  count <- charts(qr_model, k=sst, dataframe = T) %>% as.vector()
  chart <- tibble(count) %>% mutate(species = sp, age=age, sst=sst)
  return (chart)
}

## model each species in each age
## and save predicted values
data_list = list()
sp_list <- all_sp$species %>% unique()
age_list <- all_sp$age %>% unique()
i <- 1
for (sp_i in sp_list) {
  for (age_j in age_list) {

    tmp_raw <- all_sp %>% filter(age == age_j, species==sp_i) 
    tmp_predict <- predict_data(tmp_raw)
    data_list[[i]]<-tmp_predict
    i <- i + 1
    remove(tmp_raw)
    remove(tmp_predict)
  }
}

quantreg_pred_obs <- do.call("rbind", data_list)
fig3 <- quantreg_pred_obs %>% ggplot()+ geom_line(aes(x=sst, y=count, color=age),linewidth=1) +
  facet_wrap(~species,scales = "free_y") + 
  scale_color_manual(values=c("#0072B5FF", "#DC143C"),
                     labels=c("LGM","Holocene")) + 
  labs(x="Sea surface temperature (°C)", y="Abundance (#)")
fig3 <-fig3 + theme_bw() + theme(strip.background = element_blank(), strip.text = element_text(face = "italic"), legend.position = "none")

fig3%>%ggsave(file = "output/fig3.jpg", dpi = 400,width = 10, height = 8)  

fig_demonstration <- all_sp %>% filter(species=="G. ruber ruber", age=="LGM") %>% 
  ggplot()+ geom_point(aes(x=SST, y=Abundance)) + xlab("Sea Surface Temperature (C)")
fig_demonstration%>%ggsave(file = "output/ruber_points.jpg", dpi = 400,width = 3, height = 3)  
