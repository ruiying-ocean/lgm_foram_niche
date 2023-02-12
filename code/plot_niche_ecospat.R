library(ade4)
library(tidyverse)
library(ecospat)

modern <- read_csv("data/modern_4pca.csv") %>% select(-1) %>% mutate(age  = "modern") %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
lgm <- read_csv("data/lgm_4pca.csv") %>% select(-1)  %>% mutate(age  = "LGM")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
future <- read_csv("data/future4_4pca.csv") %>% select(-1) %>% mutate(age  = "future")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))

all <- rbind(modern, lgm, future) %>% pivot_longer(cols=bn:ss, names_to = "foram", values_to="biomass")

# perform PCA before scaling each column
pca.env <- all %>% select(sst:stratification) %>%
  prcomp(center=T, scale=T)

# print loading
print(pca.env$rotation)

# plot loading
# library(explor)
# explor(pca.env)
loading_plot <- factoextra::fviz_pca_var(pca.env)
ggsave("output/PCA_loading.png",loading_plot, dpi=400)

plot_niche_overlap <- function(pca, age1, age2, foram_sp,
                               title="", 
                               col1="blue", col2="red", col3="purple", 
                               tag="a", foram_threshold=1E-9){
  ## scores, check the table in below link to know loading/scores output
  ## https://aedin.github.io/PCAworkshop/articles/b_PCA.html#difference-between-covariance-based-and-correlation-based-pca-1
  scores.all.ages <- pca$x[, 1:2]
  colnames(scores.all.ages) <- c("PC1", "PC2")

  ## combine the score data with the original table
  comb.data <- cbind(all, scores.all.ages) %>% tibble()

  ## subset scores
  scores.sp.age1 <- comb.data %>% filter(age==age1 & foram==foram_sp) %>%
    filter(biomass > foram_threshold) %>% select(c(PC1, PC2))
  
  scores.sp.age2 <- comb.data %>% filter(age==age2 & foram==foram_sp) %>% 
    filter(biomass > foram_threshold) %>% select(c(PC1, PC2))

  ## PCA scores for the whole native study area
  scores.clim.age1 <- comb.data %>% filter(age==age1) %>% select(c(PC1, PC2))
  scores.clim.age2 <- comb.data %>% filter(age==age2) %>% select(c(PC1, PC2))

  
  # create grids with occurrence densities along environmental gradients.
  extend_range <- c(0,0,0,0)
  grid.clim.age1 <- ecospat.grid.clim.dyn(glob=scores.all.ages, 
                                          glob1=scores.clim.age1,
                                          sp=scores.sp.age1,
                                          R=120,
                                          th.sp=0,
                                          extend.extent=extend_range)
  
  grid.clim.age2 <- ecospat.grid.clim.dyn(glob=scores.all.ages, 
                                          glob1=scores.clim.age2,
                                          sp=scores.sp.age2,
                                          R=120,
                                          th.sp=0,
                                          extend.extent=extend_range)   
  
  eigen_value <- pca$sdev^2
  var_prop = eigen_value/sum(eigen_value)
  var_pc1 = round(var_prop[1],3) * 100
  var_pc2 = round(var_prop[2],3) * 100
  
  ## niche stability in blue, niche expansion in red, and niche unfilling in blue
  ## shade: one climate's density
  ## stability: common niche
  ## expansion: climate 2's unique niche
  ## unfilling: climate 1's unique niche
  ecospat.plot.niche.dyn(grid.clim.age1,
                         grid.clim.age2,
                         quant=0, interest=2,
                          title= title,
                          name.axis1=paste0("PC1 (", var_pc1, "%)"),
                          name.axis2=paste0("PC2 (", var_pc2, "%)"),
                          col.stab = col3,
                          col.unf = col1, col.exp=col2,
                          colZ1 =  col1, colZ2 = col2,
                         transparency = 60)
  
  ## Schoener's D index
  overlap <- ecospat.niche.overlap(grid.clim.age1, grid.clim.age2, cor = F)$D
  overlap <- round(overlap, 2)
  print(overlap)
  overlap_label <- paste0("Niche overlap: ", overlap)
  
  text(-13.5, 4, tag, font=2, cex=1.7)
  text(-10, -6, overlap_label)
  legend(x = "bottomright", legend=c(age1, age2, "overlap"), lty = c(1, 1, 1), col = c(col1, col2, col3))
}

pdf("output/optimal_niche_overlap.pdf",width=8, height=10)
  
min_biomass=1E-4

par(mfrow = c(4, 2))
lgm_color <- "#008EA0FF"
modern_color <- "orange"
future_color <- "red3"  
overlap_lgm_modern <- "brown"
overlap_modern_future <- "red1"

plot_niche_overlap(pca.env, "LGM", "modern", 
                   foram_sp="bn",
                   title = "symbiont-barren non-spinose",
                    foram_threshold = min_biomass,
                   col1=lgm_color, col2=modern_color,
                   tag="a", col3=overlap_lgm_modern)

plot_niche_overlap(pca.env, "modern", "future", 
                   foram_sp="bn",
                   title = "symbiont-barren non-spinose",
                   col1=modern_color, col2=future_color,
                   foram_threshold = min_biomass,
                   tag="b", col3=overlap_modern_future)

plot_niche_overlap(pca.env, "LGM", "modern", 
                   foram_sp="bs",
                   title =  "symbiont-barren spinose",
                   foram_threshold = min_biomass,
                   col1=lgm_color, col2=modern_color,
                   tag="c", col3=overlap_lgm_modern)

plot_niche_overlap(pca.env, "modern", "future", 
                   foram_sp="bs",
                   title = "symbiont-barren spinose",
                   foram_threshold = min_biomass,
                   col1=modern_color, col2=future_color,
                   tag="d", col3=overlap_modern_future)

plot_niche_overlap(pca.env, "LGM", "modern", 
                   foram_sp="sn",
                   title =  "symbiont-barren non-spinose",
                   foram_threshold = min_biomass,
                   col1=lgm_color, col2=modern_color,
                   tag="e", col3=overlap_lgm_modern)

plot_niche_overlap(pca.env, "modern", "future", 
                   foram_sp="sn",
                   title = "symbiont-barren non-spinose",
                   foram_threshold = min_biomass,
                   col1=modern_color, col2=future_color,
                   tag="f", col3=overlap_modern_future)

plot_niche_overlap(pca.env, "LGM", "modern", 
                   foram_sp="ss",
                   title =  "symbiont-barren spinose",
                   foram_threshold = min_biomass,
                   col1=lgm_color, col2=modern_color,
                   tag="g", col3=overlap_lgm_modern)

plot_niche_overlap(pca.env, "modern", "future", 
                   foram_sp="ss",
                   title = "symbiont-barren spinose",
                   foram_threshold = min_biomass,
                   col1=modern_color, col2=future_color,
                   tag="h", col3=overlap_modern_future)

dev.off ()


