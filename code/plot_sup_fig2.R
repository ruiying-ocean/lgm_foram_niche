## stacked barplot 

## Contact: rui.ying@bristol.ac.uk

## Read raw data
## OBSERVATION DATA (absolute and functional group format)
dir <- "https://raw.githubusercontent.com/ruiying-ocean/lgm_foram_census/main/tidy"
#dir <- "~/Science/lgm_foram_census/tidy"
lgm_fg_r <- read_csv(paste(dir, "lgm_fg_r_wsst.csv", sep = "/"))

lgm_fg_r <- lgm_fg_r %>%
  mutate(age = "LGM") %>%
  dplyr::filter(Data_Source == "margo") %>%
  rowwise()%>%
  mutate(others = sum(`symbiont-facultative non-spinose`, `symbiont-bearing spinose`,
                      `undetermined non-spinose`, `symbiont-facultative spinose`,na.rm=T)) %>%
  select(c("Latitude", "Longitude", "age", "SST", "symbiont-barren non-spinose",
           "symbiont-barren spinose", "symbiont-obligate spinose", "others"))

lgm_fg_r <- lgm_fg_r %>% pivot_longer(cols = `symbiont-barren non-spinose`:`others`, names_to = "species", values_to = "abundance")

pi_fg_r <- read_csv(paste(dir, "forcens_fg_r_wsst.csv", sep = "/")) %>%
  mutate(age = "PI") %>%
  rowwise()%>%
  mutate(others = sum(`symbiont-facultative non-spinose`, `symbiont-bearing spinose`,
                      `undetermined non-spinose`, `symbiont-facultative spinose`,na.rm=T)) %>%
  select(c("Latitude", "Longitude", "age", "SST", "symbiont-barren non-spinose",
           "symbiont-barren spinose", "symbiont-obligate spinose", "others"))

pi_fg_r <- pi_fg_r %>% pivot_longer(cols = `symbiont-barren non-spinose`:`others`, names_to = "species", values_to = "abundance")

obs_fg_r_raw <- rbind(pi_fg_r, lgm_fg_r) %>% rename(sst = SST)

## get the mid point of each bin
sst_bin_midpoints <- obs_fg_r_raw %>%
  mutate(sst_bin = cut(sst, breaks = 15)) %>%
  group_by(sst_bin) %>%
  summarize(midpoint = mean(sst, na.rm=T)) %>%
  pull(midpoint) %>% round()

## change NaN to 30
sst_bin_midpoints[16] <- 30

stacked_barplot <- obs_fg_r_raw %>%
  mutate(sst_bin = cut(sst, breaks = 15)) %>%
  ggplot(aes(x = sst_bin, y = abundance, fill = species)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_wrap(~age, ncol=2) +
  labs(fill = "Ecogroup") + 
  theme(legend.position = 'none') +
  xlab("SST Bins") +
  ylab("Average relative abundance") +
  theme_publication() +
  scale_x_discrete(labels=sst_bin_midpoints)+
  scale_fill_brewer(palette = "Spectral", direction = -1)

summary_data <- obs_fg_r_raw %>%
  mutate(sst_bin = cut(sst, breaks = 15)) %>%
  group_by(age, sst_bin, species) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE))

difference_data <- summary_data %>%
  group_by(species, sst_bin) %>%
  reframe(abundance_difference = mean_abundance[age == "PI"] - mean_abundance[age == "LGM"])

# Plot the difference in abundance
difference_plot <- ggplot(difference_data, aes(x = sst_bin, y = abundance_difference, fill = species)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(fill = "Species") +
  xlab("SST Bins") +
  ylab("Abundance difference (PI - LGM)") +
  theme_publication() +
  scale_x_discrete(labels=sst_bin_midpoints)+
  labs(fill = "Ecogroup") +
  scale_fill_brewer(palette = "Spectral", direction = -1)

library(patchwork)
combined_plot <- stacked_barplot/ difference_plot + plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = 'bold'))

ggsave("output/compositional_change.png", combined_plot,
       width = 10, height = 8, dpi = 300)
