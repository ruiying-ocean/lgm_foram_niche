library(tidyverse)
library(factoextra)
library(patchwork)
library(geometry)

## ------- PCA --------
## PCA is not good, because for every site there is all four species with difference abundance
## another way to do is to get optimal niche (e.g., abundance > 0.8) and then do the PCA

# read in niche data: only presence grids
df_modern <- read_csv("data/modern_4pca.csv") %>% select(-1) %>% mutate(age  = "preindustrial") %>% mutate(age = factor(age, levels=c("future", "preindustrial", "LGM")))
df_lgm <- read_csv("data/lgm_4pca.csv") %>% select(-1)  %>% mutate(age  = "LGM")  %>% mutate(age = factor(age, levels=c("future", "preindustrial", "LGM")))
df_future <- read_csv("data/future4_4pca.csv") %>% select(-1) %>% mutate(age  = "future")  %>% mutate(age = factor(age, levels=c("future", "preindustrial", "LGM")))
db <- rbind(df_modern, df_lgm, df_future) %>% pivot_longer(cols=bn:ss, names_to = "foram", values_to="biomass")
remove("df_modern", "df_lgm", "df_future")

get_optimal_niche <- function(foram_group, percentile = 0.75) {
 
   ## fixed percentile + different grid number in each age -> 
   ## same number of optimal sites for each group in each age despite different geographical zones
   ## simple test:
   ## db %>% filter(foram=="bn", age=="preindustrial") %>% filter(abundance > quantile(abundance, 0.8)) %>% nrow()
  return( 
    db %>% filter(foram == foram_group) %>% group_by(age) %>% 
    ## define 80th percentile for each age
    mutate(q = quantile(biomass, percentile)) %>% 
    # select the optimal abiotic niche
    filter(biomass > q) %>% select(sst:stratification, age) %>% ungroup()
  )
}

df_bn <- get_optimal_niche("bn")
df_bs <- get_optimal_niche("bs")
df_sn <- get_optimal_niche("sn")
df_ss <- get_optimal_niche("ss")

## Doing PCA transformation for environmental data
pca_bn <- df_bn %>% select(-age) %>% princomp(., cor = TRUE, scores = TRUE)
pca_bs <- df_bs %>% select(-age) %>% princomp(., cor = TRUE, scores = TRUE)
pca_sn <- df_sn %>% select(-age) %>% princomp(., cor = TRUE, scores = TRUE)
pca_ss <- df_ss %>% select(-age) %>% princomp(., cor = TRUE, scores = TRUE)

## PCA plot wrapper
plot_pca <- function(pca_object, age_array, foram_name) {
  fviz_pca_biplot(pca_object,
                  geom.ind="point",
                  col.ind = age_array,
                  pointshape=19,
                  label = "var", 
                  addEllipses = F,
                  repel = T,
                  title= foram_name, col.var="black") + 
    scale_color_manual(name = "Age",
                       labels = c("Future", "Pre-industrial", "LGM"),
                       values= c("red","orange","lightblue")) +
    theme_bw()
}

p1 <- plot_pca(pca_bn, df_bn$age, "symbiont-barren non-spinose")
p2 <- plot_pca(pca_bs, df_bs$age, "symbiont-barren spinose")
p3 <- plot_pca(pca_sn, df_sn$age, "symbiont-facultative non-spinose")
p4 <- plot_pca(pca_ss, df_ss$age, "symbiont-obligate spinose")
p <- (p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = 'A')

png("./output/optimal_niche_PCA.jpg", res=500, width = 10, height = 8, unit="in")
print(p)
dev.off()

library(ggpubr)
## plot optimal environments using box plot
ggdensity(df_bn, x="sst", fill="age", add ="mean")
ggdensity(df_bs, x="sst", fill="age", add ="mean")
ggdensity(df_sn, x="sst", fill="age", add ="mean")
ggdensity(df_ss, x="sst", fill="age", add ="mean")

ggdensity(df_bn, x="prey", fill="age", add ="mean") + xlab("Prey carbon biomass (mmol C/m3)")
ggdensity(df_bs, x="prey", fill="age", add ="mean") + xlab("Prey carbon biomass (mmol C/m3)")
ggdensity(df_sn, x="prey", fill="age", add ="mean")+ xlab("Prey carbon biomass (mmol C/m3)")
ggdensity(df_ss, x="prey", fill="age", add ="mean")+ xlab("Prey carbon biomass (mmol C/m3)")

# can't add Fe!
x <- df_bn %>% filter(age=="future") %>% select(-c(age, fe)) %>% as.matrix() %>% convhulln(., option="FA")

## Niche space using convex hull

niche_volume <- function(df, optimal) {
  df %>% group_by(age) %>% group_split() %>% 
    map(~ as.matrix(.[ , c('sst', 'po4', 'prey', 'lat', 'upwelling', 'fe', 'stratification')])) %>%
    map(scale, center=T, scale=T) %>%
    map(~ convhulln(., option="FA")) %>% 
    map_dbl("vol")
}

niche_volume(df_bn)