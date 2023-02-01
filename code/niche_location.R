library(tidyverse)
library(factoextra)
library(patchwork)
library(geometry)
library(ggthemr)
library(hypervolume) # to build robust hypervolume niche

# add PAR data
# Biomass weighted data (unit is tricky)
# Map method is different
# every run is different
# the difference is not obvious at all

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
  data4scale <- df %>% select(!all_of(except_cols))
  rest_data <-  df %>% select(all_of(except_cols))
  
  scaled_data <- scale(data4scale, scale=T) %>% data.frame() %>% tibble()
  if (combine){
    full_scaled_data <- cbind(scaled_data, rest_data)
    return(tibble(full_scaled_data))   
  } else{
    return(scaled_data)
  }
}

df_modern <- scale_data(df_modern, c("age", "bn", "bs", "sn", "ss"), combine = T)
df_future <- scale_data(df_future, c("age", "bn", "bs", "sn", "ss"), combine = T)
df_lgm <- scale_data(df_lgm, c("age", "bn", "bs", "sn", "ss"), combine = T)

global_env <- rbind(df_modern, df_lgm, df_future) %>% select(-c(bn:ss))
db_foram <- rbind(df_modern, df_lgm, df_future) %>% pivot_longer(cols=bn:ss, names_to = "foram", values_to="biomass")

## ----------------------------------------
##              subset foram group
## ----------------------------------------

# filter the minimal record
db_foram <- db_foram %>% filter(biomass > 1E-9)

all_the_bn <- db_foram %>% filter(foram == "bn") %>% select(-c(foram))
all_the_bs <- db_foram %>% filter(foram == "bs") %>% select(-c(foram))
all_the_sn <- db_foram %>% filter(foram == "sn") %>% select(-c(foram))
all_the_ss <- db_foram %>% filter(foram == "ss") %>% select(-c(foram))

## ----------------------------------------
##              get optimal subset
## ----------------------------------------

opt_bn <- all_the_bn %>% filter(biomass >= 1E-4)
opt_bs <- all_the_bs %>% filter(biomass >= 1E-4)
opt_sn <- all_the_sn %>% filter(biomass >= 1E-4)
opt_ss <- all_the_ss %>% filter(biomass >= 1E-4)


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
                        color = "black") +
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
                             values = c("#C71000FF","orange", "#008EA0FF")) +
    theme_bw() + theme(panel.grid = element_blank(),
                       text=element_text("Fira Sans", size=16),
                       panel.border = element_rect(colour = "black",
                                                   fill=NA, 
                                                   linewidth = 1))  
  
  return(p)
}

p1 <- plot_pca(opt_bn) 
p2 <- plot_pca(opt_bs)
p3 <- plot_pca(opt_sn)
p4 <- plot_pca(opt_ss)

p <- p1 + p2 + p3+ p4 + plot_annotation(tag_levels = 'a')
ggsave("output/PCA_niche.png", dpi=400, height=10, width=12)
