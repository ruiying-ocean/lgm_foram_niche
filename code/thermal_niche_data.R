library(tidyverse)
library(quantregGrowth)
library(showtext)
library(ggrepel)
library(patchwork)

global_quantlvl = 0.95

species_abbrev <- function(full_name, sep_string = ". ") {
  name_parts <- str_split(full_name, " ")[[1]]
  genus_name <- name_parts[1]
  species_name <- name_parts[2]
  
  if (length(name_parts) > 2) {
    subspecies_name <- name_parts[3]
    genus_abbrev <- str_sub(genus_name, 1, 1)
    abbrev <- paste(genus_abbrev, species_name, sep = sep_string)
    abbrev <- paste(abbrev, subspecies_name, sep=" ")
  } else {
    genus_abbrev <- str_sub(genus_name, 1, 1)
    abbrev <- paste(genus_abbrev, species_name, sep = sep_string)
  }
  
  return(abbrev)
}


## Read model csv data
load_models <- function(directory_path) {
  file_names <- list.files(directory_path, pattern = ".*_foramecogenie.csv", full.names = TRUE)
  combined_data <- data.frame()
  
  for (file_name in file_names) {
    age <- gsub(".*/(.*?)_foramecogenie.csv", "\\1", file_name)
    data <- read_csv(file_name) %>% select(-1) %>% mutate(age = age)
    combined_data <- bind_rows(combined_data, data)
  }
  
  return(combined_data)
}

## build quantile regression model
smooth_data <- function(data, x, y, quant_level = 0.9) {
  ## retrieve the name and convert to character string 
  y <- deparse(substitute(y))
  x <- deparse(substitute(x))
  data <- data %>% drop_na(all_of(c(x,y)))
  formula <- as.formula(paste(y, "~ ps(", x, ")", sep = ""))
  model <- gcrq(formula, tau = quant_level, data = data)

  fit_y <- model$fitted.values
  fit_x <- data %>% pull(x)
  
  chart <- tibble(model_y = fit_y, model_x = fit_x)
  return(chart)
}

## as above, but for a list of variables/times
loop_smooth <- function(data, i, j,...) {
  data_list <- list()
  j <- deparse(substitute(j))
  i <- deparse(substitute(i))
  vi_list <- unique(data[[i]])
  vj_list <- unique(data[[j]])
  
  n <- 1
  for (vi in vi_list) {
    for (vj in vj_list) {
##      subdata_ij <- data %>% filter({{ j }} == vj, {{ i }} == vi)
      subdata_ij <- data %>% filter(get(j) == vj, get(i) == vi)
      subsmooth_ij <- smooth_data(subdata_ij,...)
      subsmooth_ij <- subsmooth_ij %>% mutate(!!i := vi, !!j := vj)
      data_list[[n]] <- subsmooth_ij
      n <- n + 1
    }
  }
  
  combined_data <- do.call("rbind", data_list)
  return(combined_data)
}


## optimal temperature, calculated from smoothed data
thermal_opt <- function(data, long_format=TRUE){

    ## find the highest abundnace
    data <- data %>%
        group_by(species, age) %>%
        mutate(max_y = max(model_y)) %>%
        ungroup()

    ## find the corresponding temperature
    filter_data <- data %>%
        filter(model_y==max_y) %>%        
        distinct()

    if (long_format){
      return(filter_data)
    } else{
      report_data <- filter_data %>%
        pivot_wider(id_cols = species, names_from = age, values_from = model_x) %>%
        mutate(diff = PI - LGM)
      return(report_data)
    }
}


convert_to_abundance <- function(data){
    ## carbon quota source
    ## qc <- read.table("model/muffin.CBE.GIteiiaa.BASESFeTDTL_rb_foramecogem2.1/ecogem/Plankton_params.txt",header = T)
    qc <- data.frame(species = c("bn", "bs", "sn", "ss"),
                     # volume in um3
                     volume = c(1.95e+06, 2.81e+06, 3.59e+06,3.59e+06),
                     # carbon quota in mmol C/cell
                     carbon_quota_genie = c( 4.97e-06,6.85e-06, 8.51e-06, 8.51e-06))
                 
    ## every individual -> ug -> mmol C
    qc <- qc %>% mutate(carbon_quota_michaels=volume* 10/3.75E7/12*1E-3)
    data <- data %>% left_join(qc, by = c("species" = "species"))
    
    ## ind/m3 
    data <- data %>% mutate(abundance_genie = biomass/carbon_quota_genie,
                            abundance_michaels = biomass/carbon_quota_michaels)
    return(data)
}


## Read data
## OBSERVATION DATA (absolute and functional group format)
lgm_fg_a <- read_csv("~/Science/lgm_foram_census/tidy/lgm_fg_a_wsst.csv") %>% mutate(age="LGM")
lgm_fg_a <- lgm_fg_a %>% pivot_longer(cols = `Symbiont-barren Non-Spinose`:`Symbiont-obligate Spinose`,
                                  names_to = "species",values_to = "abundance")
pi_fg_a <- read_csv("~/Science/lgm_foram_census/tidy/forcens_fg_a_wsst.csv") %>% mutate(age="PI")
pi_fg_a <- pi_fg_a %>% pivot_longer(cols = `Symbiont-obligate Spinose`:`Symbiont-barren Spinose`,
                                names_to = "species",values_to = "abundance")
obs_fg_a_raw <- rbind(pi_fg_a, lgm_fg_a) %>% rename(sst=SST)
obs_fg_a_raw <- obs_fg_a_raw %>% filter(species != "Symbiont-facultative Non-Spinose")

lgm_sp <- read_csv("~/Science/lgm_foram_census/tidy/lgm_sp_a_wsst.csv") %>% mutate(age="LGM")
lgm_sp <- lgm_sp %>% pivot_longer(cols = `O. universa`:`G. crassa`,
                                  names_to = "species",values_to = "Abundance")
pi_sp <- read_csv("~/Science/lgm_foram_census/tidy/forcens_sp_a_wsst.csv") %>% mutate(age="PI")
pi_sp <- pi_sp %>% pivot_longer(cols = `D. anfracta`:`H. digitata`,
                                names_to = "species",values_to = "Abundance")
lgm_sp_list <- lgm_sp %>% pull(species) %>% unique()
pi_sp_list <- pi_sp %>% pull(species) %>% unique()
sp_list <- intersect(lgm_sp_list, pi_sp_list)
sp_list <- sp_list[!sp_list %in% c("T. humilis", "T. iota", "H. digitata",
                                   "G. uvula", "G. theyeri", "B. pumilio",
                                   "C. nitida", "D. anfracta", "H. pelagica")]

sp_list <- c("")


subset_columns <- c("SST", "species", "Abundance", "age")
obs_sp_raw <- rbind(pi_sp[subset_columns], lgm_sp[subset_columns])
obs_sp_raw <- obs_sp_raw %>% filter(species %in% sp_list)

rm(lgm_sp, pi_sp)
rm(lgm_fg_a, pi_fg_a, model_lgm, model_modern)
rm(lgm_fg_r, pi_fg_r)

## Model output
genie_fg_raw <- load_models("data/model_drived/") %>%
    select(!sn) %>%
    pivot_longer(cols=bn:ss, names_to = "species", values_to="biomass") %>%
    convert_to_abundance()%>%
    mutate(species=recode(species,
                          "bn"="Symbiont-barren Non-Spinose",
                          "bs"="Symbiont-barren Spinose",
                          "ss"="Symbiont-obligate Spinose"))

## Smooth the data
obs_fg_a_smooth <- loop_smooth(obs_fg_a_raw, i = species, j = age, x=sst, y=abundance, quant_level=global_quantlvl)
genie_fg_smooth <- loop_smooth(genie_fg_raw, i = species, j = age, x=sst, y=abundance_michaels, quant_level=global_quantlvl)
genie_fg_smooth_chl <- loop_smooth(genie_fg_raw, i = species, j = age, x=chl_total, y=abundance_michaels, quant_level=global_quantlvl)
obs_sp_smooth <- loop_smooth(obs_sp_raw, i = species, j = age, x=SST, y=Abundance, quant_level=global_quantlvl)

### -----------------------
### basic statistics

### report total data points
## obs_fg_a_raw %>% group_by(species, age) %>% summarise(n())
## obs_sp_raw %>% group_by(species, age) %>% summarise(n())


## ------------------------
## trait does not explain the difference

sp_analaysis <- thermal_opt(obs_sp_smooth) %>% group_by(age, species) %>% summarise(T_opt =model_x) %>%
  pivot_wider(id_cols = "species",names_from = "age",values_from = "T_opt") %>% mutate(diff=LGM-PI)

trait_info <- read_csv("~/Science/lgm_foram_census/fg/foram_sp_db.csv") %>%
  mutate(sp = map_vec(Name, species_abbrev)) %>% select(sp, Symbiosis, Spinose)

sp_analaysis <- merge(sp_analaysis,trait_info, by.x="species",by.y="sp")
aov(diff ~ Symbiosis * Spinose, data=sp_analaysis) %>% summary()
