library(tidyverse)
library(factoextra)


## ------- PCA --------
## PCA is not good, because for every site there is all four species with difference abundance
## another way to do is to get optimal niche (e.g., abundance > 0.8) and then do the PCA

## df <- read_csv("data/modern_4pca.csv")
df <- read_csv("data/lgm_4pca.csv")
df_bn <- df[,-c(1,3,4,5)]
df_bs <- df[,-c(1,2,4,5)]
df_sn <- df[,-c(1,2,3,5)]
df_ss <- df[,-c(1,2,3,4)]
df_bn$ecogroup <- "bn"
df_bs$ecogroup <- "bs"
df_sn$ecogroup <- "sn"
df_ss$ecogroup <- "ss"
names(df_bn)[1] <- "abundance"
names(df_bs)[1] <- "abundance"
names(df_sn)[1] <- "abundance"
names(df_ss)[1] <- "abundance"
df_bn <- df_bn %>% filter(abundance > quantile(abundance, 0.8))
df_bs <- df_bs %>% filter(abundance > quantile(abundance, 0.8))
df_sn <- df_sn %>% filter(abundance > quantile(abundance, 0.8))
df_ss <- df_ss %>% filter(abundance > quantile(abundance, 0.8))
db <- rbind(df_bn, df_bs, df_sn, df_ss)
db_simp <- db[,-c(1, 7)]
pca <- princomp(db_simp, cor = TRUE, scores = TRUE)
groups <- as.factor(db$ecogroup)

png("~/Downloads/pca_lgm.png", res=400, width = 4, height = 4, unit="in")
print(fviz_pca_biplot(pca, col.ind = groups, label = "var"))
dev.off()
