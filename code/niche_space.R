library(tidyverse)
library(factoextra)
library(patchwork)
library(geometry)
library(ggthemr)
library(hypervolume) # to build robust hypervolume niche

# add PAR data
# Biomass weighted data (unit is tricky)

## ----------------------------------------
##              Read raw data
## ----------------------------------------

df_modern <- read_csv("data/modern_4pca.csv") %>% select(-1) %>% mutate(age  = "modern") %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
df_lgm <- read_csv("data/lgm_4pca.csv") %>% select(-1)  %>% mutate(age  = "LGM")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
df_future <- read_csv("data/future4_4pca.csv") %>% select(-1) %>% mutate(age  = "future")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))

## ----------------------------------------
##              scale data
## ----------------------------------------

scale_data <- function(df, except_cols, combine=T){
  data4scale <- df %>% select(!except_cols)
  rest_data <-  df %>% select(except_cols)
  
  scaled_data <- scale(data4scale, scale=T) %>% data.frame() %>% tibble()
  if (combine){
    full_scaled_data <- cbind(scaled_data, rest_data)
    return(tibble(full_scaled_data))   
  } else{
    return(scaled_data)
  }
}

df_modern <- scale_data(df_modern, c("lat", "upwelling", "age", "bn", "bs", "sn", "ss"), combine = T)
df_future <- scale_data(df_future, c("lat", "upwelling", "age", "bn", "bs", "sn", "ss"), combine = T)
df_lgm <- scale_data(df_lgm, c("lat", "upwelling", "age", "bn", "bs", "sn", "ss"), combine = T)

global_env <- rbind(df_modern, df_lgm, df_future) %>% select(-c(bn:ss)) %>% select(!c('upwelling', 'lat'))
db_foram <- rbind(df_modern, df_lgm, df_future) %>% pivot_longer(cols=bn:ss, names_to = "foram", values_to="biomass")

## ----------------------------------------
##              subset data
## ----------------------------------------

# filter the minimal record
db_foram <- db_foram %>% filter(biomass > 1E-9)

# delete lat (duplicated info) and upwelling (binary)
db_foram <- db_foram %>% select(!c('upwelling', 'lat'))

all_the_bn <- db_foram %>% filter(foram == "bn") %>% select(-c(foram))
all_the_bs <- db_foram %>% filter(foram == "bs") %>% select(-c(foram))
all_the_sn <- db_foram %>% filter(foram == "sn") %>% select(-c(foram))
all_the_ss <- db_foram %>% filter(foram == "ss") %>% select(-c(foram))

## ----------------------------------------
##              get optimal subset
## ----------------------------------------

# get_opt_niche <- function(data, foram_group, percentile = 0.85) {
#   return( 
#     data %>% filter(foram == foram_group) %>% group_by(age) %>% 
#       ## define fixed percentile for each age
#       mutate(q = quantile(biomass, percentile)) %>% 
#       # select the optimal abiotic niche
#       filter(biomass > q) %>% ungroup()
#   )
# }

## ----------------------------------------
##              PCA plot
## ---------------------------------------- 

plot_pca <- function(data){
  # run PCA
  pca <- data %>% select(sst:stratification) %>% prcomp(scale = T, center = T)
  
  # add variables (i.e., loadings)
  PCAloadings <- data.frame(env_var = rownames(pca$rotation), pca$rotation)
  print(PCAloadings)
  
  # add group variables
  PCAvalues <- cbind(pca$x,  data[, c("biomass", "age")])
  
  # plot individual point
  p <- ggplot(PCAvalues, aes(x = PC1, y = PC2, colour = age)) +
    geom_point(aes(size=biomass))
  
  # plot loadings and annotation
  p <- p + geom_segment(data = PCAloadings, 
                        aes(x = 0, y = 0,
                            xend = (PC1*5), yend = (PC2*5)),
                        arrow = arrow(length = unit(1/2, "picas")),
                        color = "grey") +
    annotate("text",  x = (PCAloadings$PC1*5),
             y = (PCAloadings$PC2*5),
             label = PCAloadings$env_var)

  # plot variance percentage
  variance_list <- as_tibble(get_eig(pca)) %>% pull("variance.percent") %>% round(., 1)
  
  p <- p + xlab(paste0("PC1", " (", variance_list[1], "%)"))  +
    ylab(paste0("PC2", " (", variance_list[2], "%)"))
  
  # change theme
  p <- p +scale_color_manual(name = "Age",
                             labels = c("Future", "Modern", "LGM"),
                             values = c("#C71000FF","#3D3B25FF", "#008EA0FF")) +
    theme_bw() + theme(panel.grid = element_blank(),
                       text=element_text("Fira Sans", size=16),
                       panel.border = element_rect(colour = "black",
                                                   fill=NA, 
                                                   linewidth = 1))  
  
  return(p)
}

p1 <- plot_pca(all_the_bn) 
p2 <- plot_pca(all_the_bs)
p3 <- plot_pca(all_the_sn)
p4 <- plot_pca(all_the_ss)

## ----------------------------------------
##           Quantify niche space
## ---------------------------------------- 


## Niche space using convex hull, select directly related variables only
convex_hull_niche <- function(df) {
  vol <- df %>% group_by(age) %>% group_split() %>% 
    map(~as.matrix(.[ , c('sst', 'po4', 'prey', 'fe', 'mld')])) %>%
    map(scale, center=T, scale=T) %>%
    map(~ convhulln(., option="FA")) %>% 
    map_dbl("vol")
  age_order <- attributes(df$age)$levels
  return(tibble(age=age_order, niche_vol=vol))
}

## Niche space using hypervolume package, select directly related variables only
hypervolume_niche <- function(df, method="gaussian") {
  vol <- df %>% group_by(age) %>% group_split() %>% 
    map(~ as.matrix(.[ , c('sst', 'prey', 'po4', 'fe', 'mld')])) %>%
    map(~ hypervolume(., method=method)) %>% 
    map(~ get_volume(.))

  vol <- unlist(vol)
  age_order <- attributes(df$age)$levels
  return(tibble(age=age_order, niche_vol=vol))
}

# Level 1: Global climate niche
global_niche <- convex_hull_niche(global_env) %>% rename(global_niche = niche_vol)

# Level 2: All ecological niche of foram
total_niche_bn <- hypervolume_niche(all_the_bn) %>% mutate(foram="bn") %>% rename(total_niche = niche_vol)
total_niche_bs <- hypervolume_niche(all_the_bs) %>% mutate(foram="bs") %>% rename(total_niche = niche_vol)
total_niche_sn <- hypervolume_niche(all_the_sn) %>% mutate(foram="sn") %>% rename(total_niche = niche_vol)
total_niche_ss <- hypervolume_niche(all_the_ss) %>% mutate(foram="ss") %>% rename(total_niche = niche_vol)

total_niche_bn <- hypervolume_niche(all_the_bn) %>% mutate(foram="bn") %>% rename(total_niche = niche_vol)
total_niche_bs <- hypervolume_niche(all_the_bs) %>% mutate(foram="bs") %>% rename(total_niche = niche_vol)
total_niche_sn <- hypervolume_niche(all_the_sn) %>% mutate(foram="sn") %>% rename(total_niche = niche_vol)
total_niche_ss <- hypervolume_niche(all_the_ss) %>% mutate(foram="ss") %>% rename(total_niche = niche_vol)


#  Level 3: The optimal ecological niche of foram
opt_niche_bn <- hypervolume_niche(opt_bn)  %>% mutate(foram="bn")  %>% rename(optimal_niche = niche_vol)
opt_niche_bs <- hypervolume_niche(opt_bs)  %>% mutate(foram="bs") %>% rename(optimal_niche = niche_vol)
opt_niche_sn <- convex_hull_niche(opt_sn)  %>% mutate(foram="sn") %>% rename(optimal_niche = niche_vol)
opt_niche_ss <- hypervolume_niche(opt_ss) %>% mutate(foram="ss") %>% rename(optimal_niche = niche_vol)

# Combine the data into a clean dataframe
total_niche_foram <- rbind(total_niche_bn, total_niche_bs, total_niche_sn, total_niche_ss)
optimal_niche_foram <- rbind(opt_niche_bn, opt_niche_bs, opt_niche_sn, opt_niche_ss)
merged_niche_foram <- merge(total_niche_foram, optimal_niche_foram)

niche_change <- map_df(seq_len(3), ~global_niche) %>% arrange(age) %>% full_join(., merged_niche_foram) %>%
  distinct(age, foram, .keep_all=TRUE) %>%  mutate(opt_percent=optimal_niche/global_niche,
                                                   total_percent=total_niche/global_niche,)

niche_change %>% ggplot(aes(x=factor(age, levels=c("LGM","modern","future")), group=foram)) +
  geom_line(aes(y=opt_percent, color=foram))

# export to latex
# xtable::xtable(niche_change)

## ----------------------------------------
##           Plot niche change
## ---------------------------------------- 

p_niche <- niche_change %>% 
  mutate(foram = case_when(foram == 'bn' ~ 'symbiont-barren non-spinose',
                           foram == "bs" ~ 'symbiont-barren spinose',
                           foram == 'sn' ~ 'symbiont-facultative spinose',
                           foram == 'ss' ~ 'symbiont-obligate spinose')) %>%
  ggplot(aes(
  x=factor(age, levels=c("LGM", "modern", "future")),
  y=total_niche,
  group=foram,
  alpha=foram)) + 
  geom_point(aes(color=foram), size=2) + 
  geom_line(aes(color=foram), linewidth=1) + 
  xlab("") + ylab("optimal niche") +
  labs(color="") + theme_classic()+
  theme(legend.position = "none", 
        text=element_text("Fira Sans", size=8),
        panel.border = element_blank())+
  ggsci::scale_color_jco()


## ----------------------------------------
##      Incorporate niche and PCA plot
## ---------------------------------------- 


# order for inset plot location: left,bottom,right,top,
inset_width <- 0.42
inset_height <- 0.33

p_inset_bn <- p_niche +  scale_alpha_manual(values = c(1, 0.1, 0.1, 0.1))
p1<-p1 + inset_element(p_inset_bn, 0.02, 0.65, 0.02+inset_width, 0.65+inset_height, ignore_tag = F)

p_inset_bs <- p_niche +  scale_alpha_manual(values = c(0.1, 1, 0.1, 0.1))
p2<-p2 + inset_element(p_inset_bs, 0.02, 0.65, 0.02+inset_width, 0.65+inset_height, ignore_tag = F)

p_inset_sn <- p_niche +  scale_alpha_manual(values = c(0.1, 0.1, 1, 0.1))
p3 <- p3 + inset_element(p_inset_sn, 0.02, 0.3, 0.02+inset_width, 0.3+inset_height, ignore_tag = F)

p_inset_ss <- p_niche +  scale_alpha_manual(values = c(0.1, 0.1, 0.1, 1))
p4 <- p4 + inset_element(p_inset_ss, 0.02, 0.05, 0.02+inset_width, 0.05+inset_height, ignore_tag = F)

p <- p1 + p2 + p3+ p4 + plot_annotation(tag_levels = 'A')
ggsave("output/PCA_niche.png", dpi=400, height=10, width=12)
