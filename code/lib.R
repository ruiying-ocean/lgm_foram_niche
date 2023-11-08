## This script contains main functions to read the thermal performance data,
## build quantile regression model, and calculate the thermal optimum

library(tidyverse, warn.conflicts = FALSE)
library(quantregGrowth)

## abbreviate species genus names (e.g. Globigerinoides ruber to G. ruber)
species_abbrev <- function(full_name, sep_string = ". ") {
  name_parts <- str_split(full_name, " ")[[1]]
  genus_name <- name_parts[1]
  species_name <- name_parts[2]

  if (length(name_parts) > 2) {
    subspecies_name <- name_parts[3]
    genus_abbrev <- str_sub(genus_name, 1, 1)
    abbrev <- paste(genus_abbrev, species_name, sep = sep_string)
    abbrev <- paste(abbrev, subspecies_name, sep = " ")
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
    data <- read_csv(file_name) %>%
      select(-1) %>%
      mutate(age = age)
    combined_data <- bind_rows(combined_data, data)
  }

  return(combined_data)
}

## build quantile regression model for y ~ x
## quant_level is a vector of quantiles and calculate the mean/sd
smooth_qrg <- function(data, x, y, quant_level = seq(0.9, 0.99, 0.01)) {
  ## retrieve the name and convert to character string
  y <- deparse(substitute(y))
  x <- deparse(substitute(x))
  data <- data %>% drop_na(all_of(c(x, y)))
  formula <- as.formula(paste(y, "~ ps(", x, ")", sep = ""))
  model <- gcrq(formula, tau = quant_level, data = data)

  fit_y <- model$fitted.values
  fit_x <- data %>% pull(x)

  ## combine the fitted values and the original data
  chart <- cbind(fit_x, fit_y) %>% as_tibble()

  ## change name of fit_y
  fit_y_name <- paste0("model_y_", quant_level)
  names(chart) <- c("model_x", fit_y_name)

  ## get mean and sd for columns containing model_y
  model_y_columns <- chart %>% select(starts_with("model_y"))

  ## Calculate the rowwise mean for each "model_y" column
  model_y_mean <- model_y_columns %>% rowMeans()
  model_y_sd <- model_y_columns %>% apply(1, sd)

  ## merge back to the original data frame
  chart <- chart %>%
    select(-starts_with("model_y")) %>%
    cbind(model_y_mean, model_y_sd)

  return(chart)
}

## this is a nested loop function to smooth data by species and age
loop_smooth <- function(data, i, j, ...) {
  data_list <- list()
  j <- deparse(substitute(j))
  i <- deparse(substitute(i))
  vi_list <- unique(data[[i]])
  vj_list <- unique(data[[j]])

  n <- 1
  for (vi in vi_list) {
    for (vj in vj_list) {
      subdata_ij <- data %>% filter(get(j) == vj, get(i) == vi)
      subsmooth_ij <- smooth_qrg(data = subdata_ij, ...)
      subsmooth_ij <- subsmooth_ij %>% mutate(!!i := vi, !!j := vj)
      data_list[[n]] <- subsmooth_ij
      n <- n + 1
    }
  }

  combined_data <- do.call("rbind", data_list)
  return(combined_data)
}


## optimal temperature, calculated from smoothed data
thermal_opt <- function(data, long_format = TRUE) {
  ## find the highest abundnace
  data <- data %>%
    group_by(species, age) %>%
    mutate(ymax = max(model_y_mean)) %>%
    filter(model_y_mean >= 0.95 * ymax) %>%
    distinct()
  
  ## calculate the mean and sd of optimal temperature    
  report_data <- data%>%
    group_by(species, age) %>%
    summarise(opt_x_mean = mean(model_x),
              opt_x_sd = sd(model_x)) %>%
    ungroup()
  
  if (long_format) {
    return(report_data)
  } else {
    report_data <- report_data %>%
      pivot_wider(id_cols = species, names_from = age, values_from = opt_x_mean)
    return(report_data)
  }
}

## convert model biomass to abundance
convert_to_abundance <- function(data) {
  ## carbon quota source
  qc <- data.frame(
    species = c("bn", "bs", "sn", "ss"),
    # volume in um3
    volume = c(1.95e+06, 2.81e+06, 3.59e+06, 3.59e+06)
  )

  ## every individual -> ug -> mmol C
  qc <- qc %>% mutate(carbon_quota_michaels = volume * 10 / 3.75E7 / 12 * 1E-3)
  data <- data %>% left_join(qc, by = c("species" = "species"))

  ## ind/m3
  data <- data %>% mutate(abundance_michaels = biomass / carbon_quota_michaels)
  return(data)
}

## modified from https://rpubs.com/Koundy/71792
theme_publication <- function(base_size = 14, base_family = "helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size, base_family = base_family)
  + theme(
      plot.title = element_text(face = "bold"),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      #axis.title = element_text(face = "plain", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      ## legend.key = element_rect(colour = NA),
      ## legend.key.size = unit(0.2, "cm"),
      ## legend.margin = unit(0, "cm"),
      ## legend.title = element_text(face = "italic"),
      # plot.margin = unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face = "plain")
    ))
}

plot_tpc <- function(raw_data, smooth_data, x, y, se = TRUE, vline = TRUE, vlabel=FALSE, normalise_y = TRUE, colors, labels) {
  if (normalise_y) {
    ## Normalize y values to a range of 0 to 1 for raw data (group by species and age)
    ## use raw data as the reference
    scaling_factor <- smooth_data %>%
      group_by(species) %>%
      mutate(scaling_factor = 1 / max(model_y_mean + model_y_sd, na.rm = TRUE)) %>%
      select(species, age, scaling_factor) %>%
      distinct() %>%
      ungroup()

    raw_data <- raw_data %>%
      left_join(scaling_factor, by = c("species", "age")) %>%
      mutate(!!sym(y) := !!sym(y) * scaling_factor) %>%
      select(-scaling_factor)

    smooth_data <- smooth_data %>%
      group_by(species, age) %>%
      left_join(scaling_factor, by = c("species", "age")) %>%
      mutate(model_y_mean = model_y_mean * scaling_factor) %>%
      mutate(model_y_sd = model_y_sd * scaling_factor) %>%
      select(-scaling_factor)
  }

  ## plot raw data (dots)

  fig <- ggplot() +
    geom_point(data = raw_data, aes(x = !!sym(x), y = !!sym(y), color = age, shape = age), size = .3, alpha = 0.2)

  ## smoothed data (line)
  fig <- fig + geom_line(data = smooth_data, aes(x = model_x, y = model_y_mean, color = age), linewidth = 1.2)


  ## subplots by species
  if (normalise_y) {
    fig <- fig + facet_wrap(~species)
  } else {
    fig <- fig + facet_wrap(~species, scales = "free_y")
  }

  if (se) {
    fig <- fig + geom_ribbon(data = smooth_data, aes(x = model_x, ymin = model_y_mean - model_y_sd, ymax = model_y_mean + model_y_sd, fill = age), alpha = 0.3)
  }

  ## plot vertical line for optimal temperature
  if (vline) {
    fig <- fig + geom_vline(
      data = thermal_opt(smooth_data), aes(xintercept = opt_x_mean, color = age),
      linetype = "dashed", linewidth = 0.5
    )
  }
  if (vlabel) {
        fig <- fig + geom_label_repel(
      data = thermal_opt(smooth_data),
      aes(x = opt_x_mean, y = 0.9, fill = age, label = round(opt_x_mean)),
      color = "white", size = 4,
      label.r = 0.05, label.size = 0.1
    )
  }

  ## add labels
  ## if labels are not provided, skip
  if (!missing(labels)) {
      fig <- fig +
    scale_color_manual(values = colors, labels = c("LGM", "PI")) +
    scale_fill_manual(values = colors)
    }

  return(fig)
}

## a function to add global labels to a patchwork object
## source: https://github.com/thomasp85/patchwork/issues/43
add_global_label <- function(pwobj, Xlab = NULL, Ylab = NULL, Xgap = 0.03, Ygap = 0.03, ...) {
  ylabgrob <- patchwork::plot_spacer()
  if (!is.null(Ylab)) {
    ylabgrob <- ggplot() +
      geom_text(aes(x = .5, y = .5), label = Ylab, angle = 90, ...) +
      theme_void()
  }
  if (!is.null(Xlab)) {
    xlabgrob <- ggplot() +
      geom_text(aes(x = .5, y = .5), label = Xlab, ...) +
      theme_void()
  }
  if (!is.null(Ylab) & is.null(Xlab)) {
    return((ylabgrob + patchworkGrob(pwobj)) +
      patchwork::plot_layout(widths = 100 * c(Ygap, 1 - Ygap)))
  }
  if (is.null(Ylab) & !is.null(Xlab)) {
    return((ylabgrob + pwobj) +
      (xlabgrob) +
      patchwork::plot_layout(
        heights = 100 * c(1 - Xgap, Xgap),
        widths = c(0, 100),
        design = "
                                   AB
                                   CC
                                   "
      ))
  }
  if (!is.null(Ylab) & !is.null(Xlab)) {
    return((ylabgrob + pwobj) +
      (xlabgrob) +
      patchwork::plot_layout(
        heights = 100 * c(1 - Xgap, Xgap),
        widths = 100 * c(Ygap, 1 - Ygap),
        design = "
                                   AB
                                   CC
                                   "
      ))
  }
  return(pwobj)
}
