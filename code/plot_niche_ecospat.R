library(ade4)
library(tidyverse)
library("ecospat")
#library("explor")
#explor(pca.env)

modern <- read_csv("data/modern_4pca.csv") %>% select(-1) %>% mutate(age  = "modern") %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
lgm <- read_csv("data/lgm_4pca.csv") %>% select(-1)  %>% mutate(age  = "LGM")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
future <- read_csv("data/future4_4pca.csv") %>% select(-1) %>% mutate(age  = "future")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))

all <- rbind(modern, lgm, future) %>% pivot_longer(cols=bn:ss, names_to = "foram", values_to="biomass")

pca.env <- all %>% select(sst:stratification) %>%
  dudi.pca(scannf = F, nf=2, center = T, scale=T)

#pca.env$li %>% ggplot() + geom_density_2d_filled(aes(x = Axis1, y=Axis1))

plot_niche_overlap <- function(pca, data1, data2, sp, title, age1, age2, col1, col2, col3="brown", tag="a", min_biomass=1E-9){

  # variance explained (eigenvalue percentage)
  explained_variance <- pca.env$eig/sum(pca.env$eig)
  var_pc1 <- round(explained_variance[1] * 100, 2)
  var_pc2 <- round(explained_variance[2] * 100, 2)

  # loadings
  print(pca.env$c1)
  
  # score values for all ages
  scores.all.ages <- pca.env$li

  # score values for each age
  scores.sp.age1 <- suprow(pca.env, data1[which(data2[,sp]>min_biomass), 5:12])$li
  scores.sp.age2 <- suprow(pca.env, data2[which(data2[,sp]>min_biomass), 5:12])$li
  
  # PCA axis value for the whole climate
  scores.clim.age1 <- suprow(pca.env, data1[,5:12])$li
  scores.clim.age2 <- suprow(pca.env, data2[,5:12])$li
  
  # create grids with occurrence densities along environmental gradients.
  grid.clim.age1 <- ecospat.grid.clim.dyn(glob=scores.all.ages,
                                            glob1=scores.clim.age1,
                                            sp=scores.sp.age1, R=100)

   grid.clim.age2 <- ecospat.grid.clim.dyn(glob=scores.all.ages,
                                          glob1=scores.clim.age2,
                                          sp=scores.sp.age2, R=100)

  # niche stability in blue, niche expansion in red, and niche unfilling in blue
  # shade: one climate's density
  # stability: common niche
  # expansion: climate 2's unique niche
  # unfilling: climate 1's unique niche
   p <- ecospat.plot.niche.dyn(grid.clim.age1, grid.clim.age2, quant=0, interest=2,
                          title= title,
                          name.axis1=paste0("PC1 (", var_pc1, "%)"),
                          name.axis2=paste0("PC2 (", var_pc2, "%)"),
                          col.stab = col3,
                          col.unf = col1, col.exp=col2,
                          colZ1 =  col1, colZ2 = col2, transparency = 60)
   overlap <- ecospat.niche.overlap(grid.clim.age1, grid.clim.age2, cor = F)$D
   overlap <- round(overlap, 2)
   overlap_label <- paste0("Niche overlap: ", overlap)
   text(13, 3.5, overlap_label)
   text(-10, 3.5, tag, font=2, cex=1.7)
   legend(x = "bottomright", legend=c(age1, age2, "overlap"), lty = c(1, 1, 1), col = c(col1, col2, col3))
  
}

pdf("output/niche_overlap.pdf",  width=8, height=10)
par(mfrow = c(4, 2))
plot_niche_overlap(pca.env,lgm, modern, 1, title = "symbiont-barren non-spinose",
                   col1="#008EA0FF", col2="orange",
                   age1="LGM", age2="Modern", tag="a")

plot_niche_overlap(pca.env, modern, future, 1, title = "symbiont-barren non-spinose",
                   col1="orange", col2="#C71000FF",
                   age1="Modern", age2="Future",  tag="b", col3="red")

plot_niche_overlap(pca.env, lgm, modern, 2, title = "symbiont-barren spinose",
                   col1="#008EA0FF", col2="orange",
                   age1="LGM", age2="Modern",  tag="c")
plot_niche_overlap(pca.env, modern, future, 2, title = "symbiont-barren spinose",
                   col1="orange", col2="#C71000FF",col3="red",
                   age1="Modern", age2="Future",  tag="d")

plot_niche_overlap(pca.env, lgm, modern, 3, title = "symbiont-facultative non-spinose",
                   col1="#008EA0FF", col2="orange",
                   age1="LGM", age2="Modern",  tag="e")
plot_niche_overlap(pca.env, modern, future, 3, title = "symbiont-facultative non-spinose",
                   col1="orange", col2="#C71000FF",col3="red",
                   age1="Modern", age2="Future",  tag="f")

plot_niche_overlap(pca.env, lgm, modern, 4, title = "symbiont-obligate spinose",
                   col1="#008EA0FF", col2="orange",
                   age1="LGM", age2="Modern",  tag="g")
plot_niche_overlap(pca.env, modern, future, 4, title = "symbiont-obligate spinose",
                   col1="orange", col2="#C71000FF",col3="red",
                   age1="Modern", age2="Future",  tag="h")

dev.off ()


