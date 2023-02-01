library(ade4)
library(ecospat)
library(tidyverse)

modern <- read_csv("data/modern_4pca.csv") %>% select(-1) %>% mutate(age  = "modern") %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
lgm <- read_csv("data/lgm_4pca.csv") %>% select(-1)  %>% mutate(age  = "LGM")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
future <- read_csv("data/future4_4pca.csv") %>% select(-1) %>% mutate(age  = "future")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))

all <- rbind(modern, lgm, future) %>% pivot_longer(cols=bn:ss, names_to = "foram", values_to="biomass")

plot_niche_overlap <- function(data1, data2, sp, title, min_biomass=1E-9){
  pca.env <- rbind(data1,data2) %>% select(sst:stratification) %>% dudi.pca(., scannf = FALSE, nf = 2)
  # variance explained (eigenvalue percentage)
  print(pca.env$eig/sum( pca.env$eig))
  
  # correlation with variables
  print(pca.env$co)
  
  # PCA axis value for the whole study area
  scores.all.ages <- pca.env$li
  
  # PCA axis value for the species's distribution
  scores.sp.age1 <- suprow(pca.env, data1[which(data1[,sp]>min_biomass), 5:12])$li
  scores.sp.age2 <- suprow(pca.env, data2[which(data2[,sp]>min_biomass), 5:12])$li
  
  
  # PCA axis value for the whole climate
  scores.clim.age1 <- suprow(pca.env,data1[,5:12])$li
  scores.clim.age2 <- suprow(pca.env,data2[,5:12])$li
  
  # create grids with occurrence densities along environmental gradients.
  grid.clim.age1 <- ecospat.grid.clim.dyn(glob=scores.all.ages,
                                            glob1=scores.clim.age1,
                                            sp=scores.sp.age1, R=120,
                                            th.sp=0)
  
  grid.clim.age2 <- ecospat.grid.clim.dyn(glob=scores.all.ages,
                                         glob1=scores.clim.age2,
                                         sp=scores.sp.age2, R=120,
                                         th.sp=0)
  
  # niche stability in blue, niche expansion in red, and niche unfilling in blue
  # shade: one climate's density
  # stability: common niche
  # expansion: climate 2's unique niche
  # unfilling: climate 1's unique niche
  p <- ecospat.plot.niche.dyn(grid.clim.age1, grid.clim.age2, quant=0, interest=2,
                         title= title, name.axis1="PC1", name.axis2="PC2",
                         col.stab = "purple",
                         col.unf = "steelblue", col.exp="red3", 
                         colZ1 = "steelblue", colZ2 = "red3")
  overlap <- ecospat.niche.overlap(grid.clim.age1, grid.clim.age2, cor = TRUE)$D
  print(paste("overlap is:", round(overlap, 2)))
  return(pca.env)
}

png("output/niche_overlap.png", units="in", width=8, height=10, res=300)
par(mfrow = c(4, 2))
plot_niche_overlap(lgm, modern, 1, title = "symbiont-barren non-spinose, LGM & Modern")
plot_niche_overlap(modern, future, 1, title = "symbiont-barren non-spinose, Modern & Future")

plot_niche_overlap(lgm, modern, 2, title = "symbiont-barren spinose, LGM & Modern")
plot_niche_overlap(modern, future, 2, title = "symbiont-barren spinose, Modern & Future")

plot_niche_overlap(lgm, modern, 3, title = "symbiont-facultative non-spinose, LGM & Modern")
plot_niche_overlap(modern, future, 3, title = "symbiont-facultative non-spinose, Modern & Future")

plot_niche_overlap(lgm, modern, 4, title = "symbiont-obligate spinose, LGM & Modern")
plot_niche_overlap(modern, future, 4, title = "symbiont-obligate spinose, Modern & Future")

dev.off ()
