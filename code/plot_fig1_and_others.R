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
                  x = "sst", y = "abundance", label_topt = T, facet_scale="fixed",
                  colors = color_palette[1:2], labels = c("LGM", "PI"))

###### Fig 1b (data)
## not plotting two exceptional points in Symbiont-barren Spinose group
fig1b <- plot_tpc(
  raw_data = obs_fg_r_raw,
  smooth_data = obs_fg_r_smooth, x = "sst", y = "abundance",
  label_topt = T,facet_scale="fixed",
  colors = color_palette[1:2], labels = c("LGM", "PI")
)

fig1a <- fig1a + xlim(-2, 32) + ylim(-0.1,1.1)+
  theme(plot.tag = element_text(face = "bold")) +
  labs(y="Rel. abund. [model]")

fig1b <- fig1b + xlim(-2, 32) + ylim(-0.1,1.1)+
  theme(plot.tag = element_text(face = "bold")) +
  ylab("Rel. abund. [data]")

fig1a <- fig1a + theme_publication(base_size = 15) + theme(legend.position = "none")
fig1b <- fig1b + theme_publication(base_size = 15) + theme(legend.position = "none")

fig1a <- fig1a + theme(axis.title.x = element_blank())
fig1b <- fig1b + theme(axis.title.x = element_blank())

fig1 <- wrap_plots(fig1a, fig1b, ncol = 1) %>% add_global_label(
  Xlab = "Annual mean sea surface temperature (°C)",
  fontface = "plain",
  size = 5
  )
## 18 cm => 7.08 inch
## 13 cm => 5.11811 inch
ggsave(file = "output/fig1.svg", fig1, dpi = 300, width = 7.08, height = 5.11811)
system("inkscape output/fig1.svg --export-filename=output/fig1.pdf --export-dpi=300")
system("rm output/fig1.svg")

## Extended data figure 2: species level thermal performance curves
## exclude species with little abundance
exclude_sp <- c("D. anfracta", "G. uvula", "T. iota", "G. adamsi")
## remove species with little abundance
obs_sp_r_smooth <- obs_sp_r_smooth %>% filter(!species %in% exclude_sp)


ext_data_fig2 <- plot_tpc(NULL, obs_sp_r_smooth, x = "SST", y = "Abundance", label_topt = F,
                  colors = color_palette[1:2], labels = c("LGM", "PI"))
ext_data_fig2 <- ext_data_fig2 + labs(x = "Annual mean sea surface temperature (°C)", y = "Relative abundance")
ext_data_fig2 <- ext_data_fig2 + theme_publication(10) +
  theme(strip.text = element_text(face = "italic"), legend.position = "bottom") 
ext_data_fig2 %>% ggsave(file = "output/ext_data_fig2.jpg", dpi = 300, width = 7, height = 5)

### Fig3b, modelled thermal performance curves in the future
genie_fg_smooth$age <- factor(genie_fg_smooth$age, levels = c("lgm", "pi", "historical", "future1p5", "future2", "future3", "future4", '3xCO2'))

fig3b <- plot_tpc(raw_data = NULL,
                  smooth_data = filter(genie_fg_smooth, age != "historical" & age != '3xCO2'), 
                  x = "sst", 
                  y = "abundance", 
                  label_topt = T,
                  label_pos = c(seq(0.125,0.0,-0.025), seq(0.05,0.0,-0.01),c(seq(0.125,0.0,-0.025))),
                  colors = color_palette, 
                  labels = c("LGM", "PI", "2100 (+1°C)", "2100 (+2°C)", "2100 (+3°C)", "2100 (+4°C)"),
                  errorbar =F)

fig3b <- fig3b + labs(x = "Annual mean sea surface temperature (°C)", y = "Relative abundance")

fig3b <- fig3b + theme_publication() +
  theme(
      legend.position = "none",
  )

ggsave("output/fig3b.jpg",fig3b, width = 9, height = 2.8, dpi = 300)
 
# plot a PI, future4, 3xPI comparison
fig_s10 <- plot_tpc(raw_data = NULL,
                  smooth_data = filter(genie_fg_smooth, age == "future4" | age == "2.5xCO2"),
                  x = "sst", y = "abundance", errorbar = TRUE)

fig_s10 <- fig_s10 + labs(x = "Annual mean sea surface temperature (°C)", y = "Relative abundance")
fig_s10 <- fig_s10 + theme_publication() + ylim(0,1.1) +
  scale_fill_manual(labels = c("spin up", "transient"), values = c("pink", "steelblue"))+
  scale_color_manual(labels = c("spin up", "transient"), values = c("pink", "steelblue"))

ggsave("output/sup_fig10.jpg",fig_s10, width = 9, height = 2.8, dpi = 300)

# plot a PI (GMD vs. Nature version) comparison with observation as reference
model_comparison <- genie_fg_smooth %>% filter(age == "piold" | age == "pi")
pi_obs <- obs_fg_r_smooth %>% filter(age == "PI")
## modify source label for plotting
pi_obs$age <- "obs"
model_comparison <- model_comparison %>% mutate(age = ifelse(age == "piold", "Ying et al. (2023)", "This study"))
## combine both
all_pis <- rbind(model_comparison, pi_obs)

fig_s7 <- plot_tpc(raw_data = NULL,
                  smooth_data = all_pis,
                  x = "sst", y = "abundance",  facet_scale="fixed", 
                  label_topt = F, errorbar=F)

fig_s7 <- fig_s7 + labs(x = "Annual mean sea surface temperature (°C)", y = "Relative abundance")

fig_s7 <- fig_s7 + ylim(0,1.1) +
  scale_fill_manual(labels = c("Observation", "This study", "Ying et al. (2023)"), values =  c("grey",'blue','red'))+
  scale_color_manual(labels = c("Observation", "This study", "Ying et al. (2023)"), values =   c("grey",'blue','red'))
  

## remove legend title
fig_s7 <- fig_s7 + theme_publication()+
  theme(legend.title = element_blank())

ggsave("output/sup_fig7.jpg",fig_s7, width = 9, height = 2.8, dpi = 300)
