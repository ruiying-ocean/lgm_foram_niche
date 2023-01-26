library(tidyverse)
library(factoextra)
library(patchwork)
library(geometry)
library(ggthemr)

df_modern <- read_csv("data/modern_4pca.csv") %>% select(-1) %>% mutate(age  = "modern") %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
df_lgm <- read_csv("data/lgm_4pca.csv") %>% select(-1)  %>% mutate(age  = "LGM")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
df_future <- read_csv("data/future4_4pca.csv") %>% select(-1) %>% mutate(age  = "future")  %>% mutate(age = factor(age, levels=c("modern", "future", "LGM")))
db_global <- rbind(df_modern, df_lgm, df_future) %>% select(-c(bn:ss))
db_foram <- rbind(df_modern, df_lgm, df_future) %>% pivot_longer(cols=bn:ss, names_to = "foram", values_to="biomass")
remove("df_modern", "df_lgm", "df_future")

# minimal record
db_foram <- db_foram %>% filter(biomass > 1E-9)
  
get_optimal_niche <- function(data, foram_group, percentile = 0.85) {
  return( 
    data %>% filter(foram == foram_group) %>% group_by(age) %>% 
      ## define fixed percentile for each age
      mutate(q = quantile(biomass, percentile)) %>% 
      # select the optimal abiotic niche
      filter(biomass > q) %>% ungroup()
  )
}

df_bn <- get_optimal_niche(db_foram, "bn") %>% select(-c(foram,q, upwelling, lat))
df_bs <- get_optimal_niche(db_foram, "bs") %>% select(-c(foram,q, upwelling, lat))
df_sn <- get_optimal_niche(db_foram, "sn") %>% select(-c(foram,q, upwelling, lat))
df_ss <- get_optimal_niche(db_foram, "ss") %>% select(-c(foram,q, upwelling, lat))

plot_pca <- function(data){
  # run PCA
  pca <- data %>% select(sst:stratification) %>% prcomp(scale = TRUE)
  
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
                             values= c("red","grey","lightblue")) +
    theme_bw() + theme(panel.grid = element_blank(),
                       text=element_text("Fira Sans", size=16),
                       panel.border = element_rect(colour = "black",
                                                   fill=NA, 
                                                   linewidth = 1))  
  
  return(p)
}

p1 <- plot_pca(df_bn) 
p2 <- plot_pca(df_bs)
p3 <- plot_pca(df_sn)
p4 <- plot_pca(df_ss)


## Niche space using convex hull
niche_volume <- function(df) {
  vol <- df %>% group_by(age) %>% group_split() %>% 
    map(~ as.matrix(.[ , c('sst', 'po4', 'prey', 'fe','stratification')])) %>%
    map(scale, center=T, scale=T) %>%
    map(~ convhulln(., option="FA")) %>% 
    map_dbl("vol")
  age_order <- attributes(df$age)$levels
  return(tibble(age=age_order, niche_vol=vol))
}

niche_total <- niche_volume(db_global)
niche_bn <- niche_volume(df_bn)  %>% mutate(foram="bn")
niche_bs <- niche_volume(df_bs)  %>% mutate(foram="bs")
niche_sn <- niche_volume(df_sn)  %>% mutate(foram="sn")
niche_ss <- niche_volume(df_ss) %>% mutate(foram="ss")
niche_foram <- rbind(niche_bn, niche_bs, niche_sn, niche_ss)
niche_total <- map_dfr(seq_len(3), ~niche_total) %>% rename(niche_global = niche_vol)
niche_change <- merge(niche_foram, niche_total) %>% distinct(age, foram, .keep_all=TRUE) %>%
  mutate(percentage=niche_vol/niche_global)


p_niche <- niche_change %>% 
  mutate(foram = case_when(foram == 'bn' ~ 'symbiont-barren non-spinose',
                           foram == "bs" ~ 'symbiont-barren spinose',
                           foram == 'sn' ~ 'symbiont-facultative spinose',
                           foram == 'ss' ~ 'symbiont-obligate spinose')) %>%
  ggplot(aes(
  x=factor(age, levels=c("LGM", "modern", "future")),
  y=percentage,
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

# order for inset plot location: left,bottom,right,top,
inset_width <- 0.42
inset_height <- 0.33

p_inset_bn <- p_niche +  scale_alpha_manual(values = c(1, 0.1, 0.1, 0.1))
p1<-p1 + inset_element(p_inset_bn, 0.02, 0.65, 0.02+inset_width, 0.65+inset_height, ignore_tag = T)

p_inset_bs <- p_niche +  scale_alpha_manual(values = c(0.1, 1, 0.1, 0.1))
p2<-p2 + inset_element(p_inset_bs, 0.02, 0.65, 0.02+inset_width, 0.65+inset_height, ignore_tag = F)

p_inset_sn <- p_niche +  scale_alpha_manual(values = c(0.1, 0.1, 1, 0.1))
p3 <- p3 + inset_element(p_inset_sn, 0.02, 0.3, 0.02+inset_width, 0.3+inset_height, ignore_tag = F)

p_inset_ss <- p_niche +  scale_alpha_manual(values = c(0.1, 0.1, 0.1, 1))
p4 <- p4 + inset_element(p_inset_ss, 0.02, 0.05, 0.02+inset_width, 0.05+inset_height, ignore_tag = F)

p <- p1 + p2 + p3+ p4 + plot_annotation(tag_levels = 'A')
ggsave("output/PCA_niche.png", dpi=400, height=10, width=12)
## library(ggpubr)
## ## plot optimal environments using box plot
## ggdensity(df_bn, x="sst", fill="age", add ="mean")
## ggdensity(df_bs, x="sst", fill="age", add ="mean")
## ggdensity(df_sn, x="sst", fill="age", add ="mean")
## ggdensity(df_ss, x="sst", fill="age", add ="mean")

## ggdensity(df_bn, x="prey", fill="age", add ="mean") + xlab("Prey carbon biomass (mmol C/m3)")
## ggdensity(df_bs, x="prey", fill="age", add ="mean") + xlab("Prey carbon biomass (mmol C/m3)")
## ggdensity(df_sn, x="prey", fill="age", add ="mean")+ xlab("Prey carbon biomass (mmol C/m3)")
## ggdensity(df_ss, x="prey", fill="age", add ="mean")+ xlab("Prey carbon biomass (mmol C/m3)")
