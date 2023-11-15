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

color_palette <- c("#0C4876", "#699c79", "#420a68", "#932667", "#dd513a", "#fca50a")

##### Fig 1a (model)
fig1a <- plot_tpc(raw_data = genie_fg_raw %>% filter(age == "lgm" | age == "pi"), 
                  smooth_data = genie_fg_smooth %>% filter(age == "lgm" | age == "pi"),
                  x = "sst", y = "abundance_michaels", 
                 colors = color_palette[1:2], labels = c("LGM", "PI"))

###### Fig 1b (data)
## not plotting two exceptional points in Symbiont-barren Spinose group
fig1b <- plot_tpc(
  raw_data = obs_fg_a_raw %>% filter(species != "Symbiont-barren Spinose" | abundance < 600),
  smooth_data = obs_fg_a_smooth, x = "sst", y = "abundance", vline = T,
  colors = color_palette[1:2], labels = c("LGM", "PI")
)

fig1a <- fig1a + ggtitle("(a) ForamEcoGENIE Model") + xlim(-2, 32) + ylim(-0.1, 1) + theme(plot.tag = element_text(face = "bold"))
fig1b <- fig1b + ggtitle("(b) Fossil Observation") + xlim(-2, 32) + ylim(-0.1, 1) + theme(plot.tag = element_text(face = "bold"))

fig1a <- fig1a + theme_publication(base_size = 15) + theme(legend.position = "none")
fig1b <- fig1b + theme_publication(base_size = 15) + theme(legend.position = "none")

fig1a <- fig1a + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
fig1b <- fig1b + theme(axis.title.y = element_blank(), axis.title.x = element_blank())

fig1 <- wrap_plots(fig1a, fig1b, ncol = 1) %>% add_global_label(
  Ylab = "Normalised abundance",
  Xlab = "Annual mean sea surface temperature (°C)",
  fontface = "plain",
  size = 5
  )
fig1

ggsave(file = "output/fig1.svg", fig1, dpi = 300, width = 10, height = 6)
## convert to pdf
## system("inkscape output/fig1.svg --export-pdf=output/fig1.pdf")
## rsvg solution
system("rsvg-convert -f pdf output/fig1.svg > output/fig1.pdf")
## remove svg file
system("rm output/fig1.svg")

## Fig. S3, species level thermal performance curves
figs3 <- plot_tpc(obs_sp_raw, obs_sp_smooth, x = "SST", y = "Abundance", vline = FALSE, colors = color_palette[1:2], labels = c("LGM", "PI"))
figs3 <- figs3 + labs(x = "Annual mean sea surface temperature (°C)", y = "Normalised abundance")
figs3 <- figs3 + ylim(-0.2, 1)
figs3 <- figs3 + theme(strip.text = element_text(face = "italic"), legend.position = "none") + theme_publication(15)
figs3 %>% ggsave(file = "output/figs3.jpg", dpi = 400, width = 12, height = 8)

### Fig2b, modelled thermal performance curves in the future
genie_fg_smooth$age <- factor(genie_fg_smooth$age, levels = c("lgm", "pi", "historical", "future1p5", "future2", "future3", "future4", '3xCO2'))

## create an empty raw_data data.frame passing to plot_tpc
## abundance: 0, ages as in genie_fg_smooth, species as in genie_fg_smooth
fake_df <- data.frame(abundance = 0, sst = 0, age = genie_fg_smooth$age, species = genie_fg_smooth$species) %>% distinct()

fig2b <- plot_tpc(raw_data = filter(fake_df, age != "historical" & age != "3xCO2"), 
                  smooth_data = filter(genie_fg_smooth, age != "historical" & age != '3xCO2'), 
                  x = "sst", 
                  y = "abundance", 
                  vline = F, 
                  colors = color_palette, 
                  labels = c("LGM", "PI", "2100 (+1°C)", "2100 (+2°C)", "2100 (+3°C)", "2100 (+4°C)"),
                  se =F) 

fig2b <- fig2b + labs(x = "Annual mean sea surface temperature (°C)", y = "Normalised abundance")

fig2b <- fig2b + theme_publication() +
  theme(
      legend.position = "none",
  )

filter(genie_fg_smooth, grepl("future", age)) %>% thermal_opt() %>% 
  group_by(species)%>% summarise(future_mean = mean(opt_x_mean), future_sd= sd(opt_x_mean))

ggsave("output/fig2b.jpg", width = 9, height = 2.8, dpi = 300)

## plot a PI, future4, 3xPI comparison    
fig_n <- plot_tpc(raw_data = filter(fake_df, age == 'pi' | age == '3xCO2'), 
                  smooth_data = filter(genie_fg_smooth, age == "pi" | age == "3xCO2"),
                  x = "sst", y = "abundance", vline = FALSE, se = TRUE)

fig_n <- fig_n + labs(x = "Annual mean sea surface temperature (°C)", y = "Normalised abundance")
fig_n <- fig_n + theme_publication() 
fig_n
ggsave("output/fig_n.jpg", width = 9, height = 2.8, dpi = 300)
