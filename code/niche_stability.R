## This script is to test the stability of foraminiferal niche
## Strategy: if niche is stable, then given similar environmental information
## the biomass should be same (statistically, pair t-test is not significant)
## Additional constrain: different time, hight/low abundance

library(tidyverse)
rm(list = ls())

distance <- function(v1, v2){
  x = abs(mean((v1 - v2)/min(v1, v2)))
  if(any(x > 0.01)){
    return(99)
  } else{
  return(x)
  }
}

df_modern <- read_csv("data/modern_4pca.csv") %>% select(-1) %>% mutate(age  = "preindustrial") %>% mutate(age = factor(age, levels=c("future", "preindustrial", "LGM")))
df_lgm <- read_csv("data/lgm_4pca.csv") %>% select(-1)  %>% mutate(age  = "LGM")  %>% mutate(age = factor(age, levels=c("future", "preindustrial", "LGM")))
df_future <- read_csv("data/future4_4pca.csv") %>% select(-1) %>% mutate(age  = "future")  %>% mutate(age = factor(age, levels=c("future", "preindustrial", "LGM")))
db <- rbind(df_modern, df_lgm, df_future) %>% pivot_longer(cols=bn:ss, names_to = "foram", values_to="biomass")
remove("df_modern", "df_lgm", "df_future")

target_biomass <- function(df, foram, env_vector, grid_index, cutoff=0.1) {
  ## Step 1, input foram group, sampled grid, get sampled environmental information
  sampled_env <- df %>% slice(grid_index) %>% select(env_vector) %>% as.numeric()
  sampled_age <- df %>% slice(grid_index) %>% pull(age)
  
  ## Step 2, search for similar environments
  ## using the cosine similarity index (> a), a is a self-defined cut-off value
  
  sampled_df <- df %>% filter(age != sampled_age) %>% rowwise() %>%
    mutate(similarity = distance(sampled_env, c_across(env_vector))) %>%
    filter(similarity < cutoff  & similarity !=0) %>% 
    filter(similarity == min(similarity))

  if (length(sampled_df) < 1) {
    return(NULL)
  } else {
    target_biomass <- sampled_df %>% filter(similarity == max(similarity)) %>% pull(biomass) %>% mean()
    return(target_biomass)
  }
}

foram_group <- "bn"
env_names <- c("sst", "sal")
sub_db <- db %>% filter(foram==foram_group)

## data container
x <- c()
y <- c()
N <- nrow(sub_db)

## repeat for 200 times
for (i in seq(20)) {
  sample_index <- sample.int(N, 1)
  sampled_biomass <- sub_db %>% slice(sample_index) %>% pull(biomass)
  x <- append(x, sampled_biomass)
  y <- append(y, target_biomass(sub_db, foram_group, env_names, sample_index,  cutoff = 0.005))
}

## Step 4, paired t-test
t.test(x, y, paired = TRUE)

## Step 4, plot t-test
library(ggpubr)
d <- data.frame(sampled = x, target = y)
ggpaired(d, cond1 = "sampled", cond2 = "target", fill = "condition", palette = "jco")
