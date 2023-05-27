library(tidyverse)
library(quantregGrowth)
library(showtext)

model_modern <- read_csv("data/modern_4pca.csv") %>% select(-1) %>% mutate(age  = "modern")
model_lgm <- read_csv("data/lgm_4pca.csv") %>% select(-1)  %>% mutate(age  = "LGM")
genie_raw <- rbind(model_modern, model_lgm) %>% pivot_longer(cols=bn:ss, names_to = "species", values_to="biomass")

## observation data
lgm_fg <- read_csv("~/foram_core/tidy/lgm_fg_a_wsst.csv") %>% mutate(age="LGM")
lgm_fg <- lgm_fg %>% pivot_longer(cols = `Symbiont-barren Non-Spinose`:`Symbiont-obligate Spinose`,
                                  names_to = "species",values_to = "abundance")
pi_fg <- read_csv("~/foram_core/tidy/forcens_fg_a_wsst.csv") %>% mutate(age="PI")
pi_fg <- pi_fg %>% pivot_longer(cols = `Symbiont-obligate Spinose`:`Symbiont-barren Spinose`,
                                names_to = "species",values_to = "abundance")
obs_raw <- rbind(pi_fg, lgm_fg) %>% rename(sst=SST)

smooth_data <- function(data, x, y, quant_level = 0.9) {

  ## retrieve the name and convert to character string 
  y <- deparse(substitute(y))
  x <- deparse(substitute(x))
  data <- data %>% drop_na(all_of(c(x,y)))
  formula <- as.formula(paste(y, "~ ps(", x, ")", sep = ""))
  model <- gcrq(formula, tau = quant_level, data = data)

  fit_y <- model$fitted.values
  fit_x <- data %>% pull(x)
  
  chart <- tibble(model_y = fit_y, model_x = fit_x)
  return(chart)
}

# test <- obs_raw %>% filter(species=="Symbiont-obligate Spinose", age=="LGM")
# smooth_data(test, x=sst,y=abundance)


loop_smooth <- function(data, i,j,...) {
  data_list <- list()
  j <- deparse(substitute(j))
  i <- deparse(substitute(i))
  vi_list <- unique(data[[i]])
  vj_list <- unique(data[[j]])
  
  n <- 1
  for (vi in vi_list) {
    for (vj in vj_list) {
##      subdata_ij <- data %>% filter({{ j }} == vj, {{ i }} == vi)
      subdata_ij <- data %>% filter(get(j) == vj, get(i) == vi)
      subsmooth_ij <- smooth_data(subdata_ij,...)
      subsmooth_ij <- subsmooth_ij %>% mutate(!!i := vi, !!j := vj)
      data_list[[n]] <- subsmooth_ij
      n <- n + 1
    }
  }
  
  combined_data <- do.call("rbind", data_list)
  return(combined_data)
}


obs_smooth <- loop_smooth(obs_raw, i = species, j = age, x=sst, y=abundance, quant_level=0.9)
genie_smooth <- loop_smooth(genie_raw, i = species, j = age, x=sst, y=biomass, quant_level=0.9)


# Set the global theme for text properties
theme_set(
  theme_bw() +
    theme(
      text = element_text(
        family = "Roboto",    # Font family
        size = 15             # Font size
      )
    )
)

fig1b <- ggplot() +  
  geom_point(data=genie_raw, aes(x=sst, y =biomass, color=age), size=0.1, alpha=0.5)+
  geom_line(data=genie_smooth, aes(x=model_x, y =model_y, color=age),  linewidth=1)+
  facet_wrap(~species, scale="free_y",nrow=1)
fig1b <- fig1b +   scale_color_manual(values=c("#0072B5FF", "#DC143C"), labels=c("LGM","Holocene"))  + 
  labs(x="Sea surface temperature (°C)", y="Biomass (mmol C/m3)")
fig1b <- fig1b + theme(strip.background = element_blank(), strip.text = element_blank(), legend.position = "none")

fig1a <- ggplot()+
  geom_point(data=obs_raw,aes(x=sst, y=abundance, color=age), size=0.1, alpha=0.5) +
  geom_line(data=obs_smooth,aes(x=model_x, y=model_y, color=age),linewidth=1) + 
  facet_wrap(~species,scales = "free_y",nrow=1)

fig1a <- fig1a +  scale_color_manual(values=c("#0072B5FF", "#DC143C"),
                                 labels=c("LGM","Holocene")) + labs(x="", y="Abundance (#)")
fig1a <- fig1a + theme(strip.background = element_blank(), strip.text = element_text(face = "italic"), legend.position = "none")
fig1 <- ggpubr::ggarrange(fig1a, fig1b, nrow = 2)

fig1 %>% ggsave(., file="output/fig1.jpg", dpi=500, width=12, height = 5)

###### Fig 3
lgm_sp <- read_csv("~/foram_core/tidy/lgm_sp_a_wsst.csv") %>% mutate(age="LGM")
lgm_sp <- lgm_sp %>% pivot_longer(cols = `O. universa`:`G. theyeri`,
                                  names_to = "species",values_to = "Abundance")
pi_sp <- read_csv("~/foram_core/tidy/forcens_sp_a_wsst.csv") %>% mutate(age="PI")
pi_sp <- pi_sp %>% pivot_longer(cols = `D. anfracta`:`H. digitata`,
                                names_to = "species",values_to = "Abundance")
# top_lgm_sp <- lgm_sp %>% group_by(species) %>% summarise(ab = mean(Abundance, na.rm=T))%>%
#   arrange(desc(ab)) %>% head(30) %>% pull(species)
# top_pi_sp <- pi_sp %>% group_by(species) %>% summarise(ab = mean(Abundance, na.rm=T))%>%
#   arrange(desc(ab)) %>% head(30) %>% pull(species)
# select_species <- intersect(top_lgm_sp, top_pi_sp) 


select_species <- c("G. bulloides", "N. pachyderma",
                    "N. dutertrei", "N. incompta", "G. inflata", "G. glutinata",
                    "G. ruber albus", "G. ruber ruber","G. truncatulinoides",
                    "P. obliquiloculata", "T. quinqueloba","T. sacculifer")

subset_columns <- c("SST", "species", "Abundance", "age")
obs_sp_raw <- rbind(pi_sp[subset_columns], lgm_sp[subset_columns])
obs_sp_raw <- all_sp %>% filter(species %in% select_species)

obs_sp_smooth <- loop_smooth(obs_sp_raw, i = species, j = age, x=SST, y=Abundance, quant_level=0.9)

fig3 <- ggplot()+ 
  geom_point(data=obs_sp_raw, aes(x=SST, y =Abundance, color=age), size=0.1, alpha=0.5)+
  geom_line(data=obs_sp_smooth,aes(x=model_x, y=model_y, color=age),linewidth=1) +
  facet_wrap(~species,scales = "free_y") + 
  scale_color_manual(values=c("#0072B5FF", "#DC143C"),
                     labels=c("LGM","Holocene")) + 
  labs(x="Sea surface temperature (°C)", y="Abundance (#)")
fig3 <-fig3  + theme(strip.background = element_blank(), strip.text = element_text(face = "italic"), legend.position = "none")

fig3%>%ggsave(file = "output/fig3.jpg", dpi = 400,width = 10, height = 8)  

