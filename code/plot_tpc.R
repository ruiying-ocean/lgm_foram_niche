## This script contains functions to plot Fig. 1 and Fig. S3

## Contact: rui.ying@bristol.ac.uk

## load visualisation packages
library(showtext)
library(ggrepel)
library(patchwork)
source("code/lib.R")

## read presaved Rdata
load("data/obs_smooth.Rdata")
load("data/obs_raw.Rdata")
load("data/genie_fg_smooth.Rdata")
load("data/genie_fg_raw.Rdata")

# Set the global theme for text properties
##### Fig 1a (model)
fig1a <- plot_tpc(raw_data = genie_fg_raw %>% filter(age == "lgm" | age == "pi"), smooth_data = genie_fg_smooth %>% filter(age == "lgm" | age == "pi"), x = "sst", y = "abundance_michaels")
fig1a <- fig1a+
    labs(
        x = "Sea surface temperature (°C)",
        y = expression("Abundance (" * "#/m"^3 * ")")
    )


###### Fig 1b (data)
fig1b <- plot_tpc(raw_data = obs_fg_a_raw, smooth_data = obs_fg_a_smooth, x = "sst", y = "abundance")
fig1b <- fig1b + labs(x = "Sea surface temperature (°C)", y = "Abundance (#)")

fig1a <- fig1a + ggtitle("(a) ForamEcoGENIE Model") + xlim(-2, 32) + theme(plot.tag = element_text(face = "bold"))
fig1b <- fig1b + ggtitle("(b) Fossil Observation") + xlim(-2, 32) + theme(plot.tag = element_text(face = "bold"))

fig1a <- fig1a + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
fig1b <- fig1b + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))

fig1a <- fig1a +  theme_publication(base_size = 15) + theme(legend.position = "none")
fig1b <- fig1b +  theme_publication(base_size = 15) + theme(legend.position = "none")

fig1 <- fig1a / fig1b

fig1 %>% ggsave(., file = "output/fig1.png", dpi = 300, width = 10, height = 6)


figs3 <- plot_tpc(obs_sp_raw, obs_sp_smooth, x = "SST", y = "Abundance")
figs3 <- figs3 + theme(strip.text = element_text(face = "italic"), legend.position = "none")
figs3 %>% ggsave(file = "output/figs3.jpg", dpi = 400, width = 12, height = 8)

### Fig2b
## the same color as the python script
color_palette <- c(c_cold, c_warm, "#420a68", "#932667", "#dd513a", "#fca50a")

genie_fg_smooth$age <- factor(genie_fg_smooth$age, levels = c("lgm", "pi", "historical", "future1p5", "future2", "future3", "future4"))

fig2b <- genie_fg_smooth %>%
  filter(age != "historical") %>%
  ggplot(aes(x = model_x, y = model_y, color = age)) +
  geom_line(linewidth = 1) +
  facet_wrap(~species, scales = "free_y", nrow = 1) +
  labs(
    x = "Sea surface temperature (°C)",
    y = expression("Abundance (" * "#/m"^3 * ")"), # LaTeX expression for the y-axis label
    color = ""
  )

fig2b <- fig2b +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.text.y = element_text(size = 8), # Adjust the font size here (smaller value)
    strip.text = element_text(size = 12),
    strip.background = element_blank()
  ) +
  scale_color_manual(
    values = color_palette,
    labels = c("Last Glacial Maximum", "Pre-industrial", "2100 (+1°C)", "2100 (+2°C)", "2100 (+3°C)", "2100 (+4°C)")
  )

ggsave("output/fig2b.jpg", width = 9, height = 2.5, dpi = 300)

## as fig1 but plot chl
figs5 <- ggplot() +
  geom_line(data = genie_fg_smooth_chl %>% filter(age != "historical"), aes(x = model_x, y = model_y, color = age), linewidth = 1)
figs5 <- figs5 + facet_wrap(~species, scale = "free_y", nrow = 1)
figs5 <- figs5 +
  scale_color_manual(
    values = color_palette,
    labels = c("Last Glacial Maximum", "Pre-industrial", "2100 (+1.5°C)", "2100 (+2°C)", "2100 (+3°C)", "2100 (+4°C)")
  ) +
  labs(
    x = expression("Total Chlorophyll (" * "mg/m"^3 * ")"),
    y = expression("Abundance (" * "#/m"^3 * ")")
  )
figs5 <- figs5 + theme(legend.position = "bottom")
ggsave("output/figs5.jpg", figs5, width = 8, height = 4, dpi = 300)
