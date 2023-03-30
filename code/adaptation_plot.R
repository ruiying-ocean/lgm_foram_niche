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
        summarise(mean_count = mean(Count, na.rm = T))
    
    return (data)
}

thermal_trait <- function(data){
    data <- data %>% filter(mean_count >= 1)
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
        filter(mean_count == max(mean_count, na.rm = T)) %>%
        distinct() %>%
        rename(ct_opt = bin_mid_pnt) %>% select(!mean_count)

                                        # merge two data frames
    return(merge(ct_max_min, ct_opt))
}



bin_niche <- all_niche %>% bin_along_var(., "SST", n=45)

p1 <- bin_niche %>% ggplot(aes(x=bin_mid_pnt, y=mean_count)) + geom_point(aes(color=age),size=1)+
    geom_line(aes(color=age),linewidth=.85, alpha=0.8)+facet_wrap(~Species, scale="free_y") 
p1 <- p1 +  scale_color_manual(values=c("#0072B5FF", "#E18727FF", "#7876B1FF"),
                     labels=c("LGM","PI", "Modern")) + labs(x="Sea surface temperature", y="Abundance (#)")

trait_df <- thermal_trait(bin_niche)
trait_df <- trait_df %>%  pivot_longer(cols = ct_min:ct_opt, names_to = 'trait')

trait_diff <- trait_df %>% group_by(trait, age) %>% summarise(mean=mean(value)) %>% 
  pivot_wider(names_from = c("age"), values_from = mean) %>% 
  mutate(diff = pi-lgm)

trait_df %>% ggplot(aes(x=trait, y=value)) + geom_boxplot(aes(fill=age))

ggboxplot(trait_df, x="age", y="value", fill="age",palette = "jco", facet.by = "trait") +
  stat_compare_means(paired=T)


## Species order as functional groups
## symbiont-barren non-spinose
## g1 <- c("N. pachyderma", "N. incompta", "T. quinqueloba")
## symbiont-barren spinose
## g2 <- c("G. bulloides")
## symbiont-falcultative
## g3 <- c("G. inflata", "N. dutertrei", "G. glutinata","G. menardii")
## symbiont-obligate
## g4 <- c("G. ruber p", "G. ruber w", "T. sacculifer", "O. universa")

## tmp <- tmp %>% mutate( function_group = case_when(
##   Species %in% g1 ~ "symbiont-barren non-spinose",
##   Species %in% g2 ~ "symbiont-barren spinose",
##   Species %in% g3 ~ "symbiont-falcultative non-spinose",
##   Species %in% g4 ~ "symbiont-obligate spinose",
## ))
