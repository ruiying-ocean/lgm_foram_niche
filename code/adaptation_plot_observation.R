library(santoku)
library(patchwork)
library(ggpubr)
source("code/adaptation_data.R")

bin_along_var <- function(data, var, n, ...){
    env_var = data %>% pull(var)
    data <- data %>% mutate(bin_lvl=chop_evenly(env_var, n))
    step <- (max(env_var) - min(env_var))/(n)
    bin_mid_pnt <- min(env_var) + step/2 * seq(1, 2*n, 2)
    bin_lvl <- tab_evenly(env_var, n) %>% names()
    bin_df <- data.frame(bin_mid_pnt, bin_lvl)
    data <- data %>% left_join(bin_df)
    
    ## take average for abundance in each bin
    data <- data %>% group_by(Species, age, bin_mid_pnt) %>%
        summarise(mean_abundance = mean(Count, na.rm = T))
    
    return (data)
}

thermal_trait <- function(data){
    data <- data %>% filter(mean_abundance >= 1)
    ## calculate Tmin, Tmax
    ct_max_min <- data  %>%
        group_by(Species, age) %>%
        mutate(ct_min = min(bin_mid_pnt, na.rm = T),
               ct_max = max(bin_mid_pnt, na.rm = T)) %>% 
        select(c(Species, age, ct_min, ct_max)) %>%
        distinct() %>%
        ungroup()
    
    ## calculate Topt
    ct_opt <- data %>% group_by(Species, age) %>% 
        filter(mean_abundance == max(mean_abundance, na.rm = T)) %>%
        distinct() %>%
        rename(ct_opt = bin_mid_pnt) %>% select(!mean_abundance)

    ## merge two data frames
    return(merge(ct_max_min, ct_opt))
}



bin_niche <- niche_data %>% bin_along_var(., "SST", n=45)

trait_df <- thermal_trait(bin_niche)

trait_df <- trait_df %>% mutate(width=ct_max-ct_min)
# Species order as functional groups
# symbiont-barren non-spinose
g1 <- c("N. pachyderma", "N. incompta", "T. quinqueloba")
# symbiont-barren spinose
g2 <- c("G. bulloides")
# symbiont-falcultative
g3 <- c("G. inflata", "N. dutertrei", "G. glutinata","G. menardii")
# symbiont-obligate
g4 <- c("G. ruber p", "G. ruber w", "T. sacculifer", "O. universa")

trait_df <- trait_df %>% mutate(function_group = case_when(
  Species %in% g1 ~ "symbiont-barren non-spinose",
  Species %in% g2 ~ "symbiont-barren spinose",
  Species %in% g3 ~ "symbiont-falcultative non-spinose",
  Species %in% g4 ~ "symbiont-obligate spinose",
), .after = age)

trait_dfl <- trait_df %>%  pivot_longer(cols = ct_min:width, names_to = 'trait')

# functional group difference in trait change
trait_df %>% group_by(Species)%>% mutate(diff = ct_max-lag(ct_max)) %>% ungroup() %>%
  group_by(function_group) %>% summarise(mean(diff, na.rm = T))

trait_df %>% group_by(Species)%>% mutate(diff = ct_opt-lag(ct_opt)) %>% ungroup() %>%
  group_by(function_group) %>% summarise(mean(diff, na.rm = T))

trait_mean_diff <- trait_dfl %>% group_by(age, trait) %>% summarise(mean=mean(value)) %>% 
  pivot_wider(names_from = age, values_from = mean)  %>% mutate(diff=pi-lgm)

trait_sp_diff <- trait_dfl %>% group_by(Species,age, trait) %>% summarise(mean=mean(value)) %>% 
  pivot_wider(names_from = age, values_from = mean)  %>% mutate(diff=pi-lgm)

p1 <- bin_niche %>% ggplot(aes(x=bin_mid_pnt, y=mean_abundance)) + geom_point(aes(color=age),size=1)+
    geom_line(aes(color=age),linewidth=.85, alpha=0.8)+facet_wrap(~Species, scale="free_y") 

p1 <- p1 +  scale_color_manual(values=c("#0072B5FF", "#E18727FF", "#7876B1FF"),
                     labels=c("LGM","PI", "Modern")) + labs(x="Sea surface temperature (°C)", y="Abundance (#)")

p1 <- p1+ theme_bw() + theme(strip.background = element_blank(), strip.text = element_text(face = "italic"))

p2 <- ggboxplot(trait_dfl, x="age", y="value", fill="age", facet.by = "trait",
                palette = c("#0072B5FF", "#E18727FF", "#7876B1FF"),
                width=0.4,
                orientation="horizontal",  nrow=1,
                panel.labs = list(trait=c("Thermal maximum","Thermal minimum","Thermal optimum", "Thermal width")))

p2 <- p2 +  stat_compare_means(paired=T, label = "p.format", 
                               label.x = 0.5, label.y = 22,
                               size = 3)

p2 <- p2 + labs(x="", y="Sea surface temperature (°C)") + guides(fill="none")
p2 <- p2+theme_bw()
p <- p1/p2+plot_layout(heights = c(5, 1))

ggsave("output/Niche_SST_adaptation.png", p, width = 12, height = 10, dpi=400)


