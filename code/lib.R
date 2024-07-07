## This script contains main functions to read the thermal performance data,
## build quantile regression model, and calculate the thermal optimum

## contact: rui.ying@bristol.ac.uk

library(tidyverse, warn.conflicts = FALSE)
library(quantregGrowth)
library(ggpubr)

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

## this is a nested loop function to run smooth_qrg for various species and ages
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


## report the thermal optimum min, max, mean, sd
## it is defined as the temperature at which the biomass is 50% of the maximum (adjustable in `Topt_coef`)
## It keep the same standard in different ages
thermal_opt <- function(data, Topt_coef = 0.5, long_format = TRUE) {
    ymax <- data %>%
    group_by(species) %>%
    slice(which.max(model_y_mean)) %>%
    select(species, model_y_mean) %>%
    rename(ymax = model_y_mean)

  data <- left_join(data, ymax, by = c("species"))

  data <- data %>%
    filter(model_y_mean >= Topt_coef * ymax)

  ## calculate the mean and sd of optimal temperature
  report_data <- data %>%
    group_by(species, age) %>%
    summarise(
      Topt_mean = mean(model_x, na.rm = TRUE),
      Topt_sd = sd(model_x, na.rm = TRUE),
      Topt_min = min(model_x, na.rm = TRUE),
      Topt_max = max(model_x, na.rm = TRUE)
    ) %>%
    ungroup()

  if (long_format) {
    return(report_data)
  } else {
    report_data <- report_data %>% pivot_wider(
      names_from = age,
      values_from = c("Topt_mean", "Topt_sd", "Topt_min", "Topt_max"),
      names_glue = "{age}_{.value}"
    )
    return(report_data)
  }
}

## modified from https://rpubs.com/Koundy/71792
theme_publication <- function(base_size = 14, base_family = "helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size, base_family = base_family)
  + theme(
      plot.title = element_text(face = "plain"),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      # axis.title = element_text(face = "plain", size = rel(1)),
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

plot_tpc <- function(raw_data, smooth_data, x, y, errorbar = TRUE, label_topt = TRUE, facet_scale="free_y",label_pos, colors, labels, linetype) {
  fig <- ggplot()

  ## plot raw data (dots)
  if (!is.null(raw_data)) {
    fig <- fig +
      geom_point(data = raw_data, aes(x = !!sym(x), y = !!sym(y), color = age, shape = age), size = .2, alpha = 0.15)
  }

  ## smoothed data (line)  
  ## if specify the linetype
  if (!missing(linetype)) {
        fig <- fig + geom_line(data = smooth_data, aes(x = model_x, y = model_y_mean, color = age, linetype = age), linewidth = 0.7)
  } else{
      ## by default, use solid line
      fig <- fig + geom_line(data = smooth_data, aes(x = model_x, y = model_y_mean, color = age), linewidth = 0.7)
  }

  ## subplot by species
  fig <- fig + facet_wrap(~species, scales = facet_scale)

  ## plot the ensemble standard deviation
  if (errorbar) {
    fig <- fig + geom_ribbon(data = smooth_data, aes(x = model_x, ymin = model_y_mean - model_y_sd, ymax = model_y_mean + model_y_sd, fill = age), alpha = 0.2)
  }

  ## plot vertical line for mean optimal temperature
  if (label_topt) {
    opt_data <- thermal_opt(smooth_data)
    print(opt_data)

    if (missing(label_pos)) {
      ## default for the LGM and PI main figure
      opt_data <- opt_data %>% mutate(y.pos = case_when(
        toupper(age) == "LGM" ~ -0.1,
        toupper(age) == "PI" ~ -0.05,
      ))
    } else {
        opt_data <- opt_data %>% arrange(species, age)
        opt_data <- opt_data %>% mutate(y.pos = label_pos)
    }

    fig <- fig +
      geom_point(
        data = opt_data,
        aes(
          x = Topt_mean,
          y = y.pos,
          color = age,
        )
      ) +
      geom_segment(
        data = opt_data,
        aes(
          x = Topt_min, xend = Topt_max,
          y = y.pos, yend = y.pos,
          group = age, color = age
        ),
        linewidth = 0.5,
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


## helper function to return lm coefficients as a list
lm_coeffs <- function(x, y) {
  coeffs = as.list(coefficients(lm(y~x)))
  names(coeffs) = c('i', "s")
  return(coeffs)
}

boot_lm <- function(data, x, y){
  library(data.table)
  set.seed(999999)
  ## generate bootstrap samples of slope ('s') and intercept ('i')
  nboot <- 1000
  mtboot <- lapply(seq_len(nboot), function(i) {
    ## resample
    sampled_data <- data[sample(nrow(data), replace = TRUE), ]
    ## get coefficients of linear model
    lm_coeffs(sampled_data[[x]], sampled_data[[y]])
  })
  mtboot <- rbindlist(mtboot)
  return(mtboot)
}

plot_lm <- function(data, x, y, label_pos,...){
  # Use R2 instead of R
  p <- ggscatter(data, x = x, y = y, add = "reg.line",...) +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
             label.x = label_pos[1], label.y = label_pos[2], size = 4)

  ## plot the bootstrapped regression results
  mtboot <- boot_lm(data, x, y)
  p <- p + geom_abline(aes(intercept=i, slope=s), data = mtboot, linewidth=0.1, color='grey', alpha=0.1)
  p <- p + geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2)

  ## annotate slope
  ## 95% CI
  ci_lower <- quantile(mtboot$s, 0.025) %>% round(1)
  ci_upper <- quantile(mtboot$s, 0.975) %>% round(1)
  median_value  <- quantile(mtboot$s, 0.5) %>% round(1)
  p <- p + annotate("text", x = label_pos[3], y = label_pos[4],  size = 4,
                    label = paste0("slope: ",median_value, " (", ci_lower,", ", ci_upper,")"))
  return(p)
}
