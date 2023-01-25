library(tidyverse)
library(factoextra)
library(patchwork)
library(geometry)

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

df_bn <- get_optimal_niche(db_foram, "bn") %>% select(-c(foram,q, upwelling))
df_bs <- get_optimal_niche(db_foram, "bs") %>% select(-c(foram,q, upwelling))
df_sn <- get_optimal_niche(db_foram, "sn") %>% select(-c(foram,q, upwelling))
df_ss <- get_optimal_niche(db_foram, "ss") %>% select(-c(foram,q, upwelling))

plot_pca <- function(data){
  pca <- data %>% select(sst:stratification) %>% prcomp(scale = TRUE)
  basic_plot <- fviz_pca_ind(pca, label="none")
  advanced_data <- cbind(basic_plot$data, data[, c("biomass", "age")])
  p <- ggplot(advanced_data, aes(x=x,y=y,col=age)) + geom_point(aes(size=biomass)) +
    scale_color_manual(name = "Age",
                       labels = c("Future", "Pre-industrial", "LGM"),
                       values= c("red","orange","lightblue")) +  theme_bw()
  return(p)
}

plot_pca(df_bn)
plot_pca(df_bs)
plot_pca(df_sn)
plot_pca(df_ss)

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

library(ggthemr)

p_niche <- niche_change %>% 
  mutate(foram = case_when(foram == 'bn' ~ 'symbiont-barren non-spinose',
                           foram == "bs" ~ 'symbiont-barren spinose',
                           foram == 'sn' ~ 'symbiont-facultative spinose',
                           foram == 'ss' ~ 'symbiont-obligate spinose')) %>%
  ggplot(aes(
  x=factor(age, levels=c("LGM", "modern", "future")),
  y=percentage,
  group=foram)) + 
  geom_point(aes(color=foram), size=4) + 
  geom_line(aes(color=foram), linewidth=1.2) + 
  xlab("") + ylab("Optimal niche percentage") +
  labs(color="")+
  theme(legend.position = "top", text=element_text("Fira Sans", size=18),
        panel.border = element_rect(fill = NA, linewidth=1),
        panel.background = element_rect(fill = NA, color = "#F9F6EE"),
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'dashed',
                                        colour = "grey"))+
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggsci::scale_color_jco()

p_niche
ggsave("output/niche_change.png",p_niche, dpi=300,width = 8, height = 6)
