source("code/thermal_niche_data.R")
### Plot !

# Set the global theme for text properties
theme_set(
  theme_bw() +
    theme(
      text = element_text(
        size = 15             # Font size
      )
    )
)

c_cold <- "#0C4876"
c_warm <- "#98bad5"

###### Fig 1
## plot raw data (dots)
fig1a <- ggplot()+ geom_point(data=obs_fg_a_raw, aes(x=sst, y=abundance, color=age), size=1.5, alpha=0.1)
## smoothed data (line)
fig1a <- fig1a + geom_line(data=obs_fg_a_smooth,aes(x=model_x, y=model_y, color=age),linewidth=1)
## subplots by species
fig1a <- fig1a + facet_wrap(~species,scales = "free_y",nrow=1)
fig1a <- fig1a + theme(legend.position = 'none', legend.background=element_blank())
## plot vertical line for optimal temperature
fig1a <- fig1a + geom_vline(data = thermal_opt(obs_fg_a_smooth), aes(xintercept = model_x, color = age), linetype = "dashed", linewidth=0.5)

## change color and labels
fig1a <- fig1a +  scale_color_manual(values=c(c_cold,c_warm), labels=c("LGM","PI")) + labs(x="Sea surface temperature (°C)", y="Abundance (#)")
## miscelaneous theme settings

fig1a<-fig1a + geom_label_repel(data = thermal_opt(obs_fg_a_smooth),
                                aes(x = model_x, y=model_y, fill = age, label=round(model_x)),
                                color="white", size=4, 
                                nudge_y = 30, label.r = 0.05, label.size=0.1)+
  scale_fill_manual(values=c(c_cold,c_warm))
fig1a
## do the same for GENIE model output
## subset the LGM and PI data
fig1b <- ggplot() + geom_point(data=genie_fg_raw %>% filter(age=="lgm" | age=="pi"), aes(x=sst, y =abundance_michaels, color=age), size=1, alpha=0.2)
fig1b <- fig1b + geom_line(data=genie_fg_smooth %>% filter(age=="lgm" | age=="pi"), aes(x=model_x, y =model_y, color=age),  linewidth=1)
fig1b <- fig1b + facet_wrap(~species, scale="free_y",nrow=1)

fig1b <- fig1b  + scale_color_manual(values=c(c_cold,c_warm), labels=c("LGM","Pre-industrial")) +
  labs(x="Sea surface temperature (°C)",
       y = expression("Abundance (" * "#/m"^3 * ")"))

fig1b <- fig1b + theme(legend.position = "none")
fig1b <- fig1b + geom_vline(data = thermal_opt(genie_fg_smooth %>% filter(age=="lgm" | age=="pi")), aes(xintercept = model_x, color = age), linetype = "dashed", linewidth=0.5)
fig1b <- fig1b + scale_y_continuous(limits = ~ c(min(.x), ceiling(max(.x)))*1.5)
fig1b <- fig1b + geom_label_repel(data = thermal_opt(genie_fg_smooth %>% filter(age=="lgm" | age=="pi")),
                                aes(x = model_x, y=model_y, fill = age, label=round(model_x)),
                                color="white",
                                size=4,
                                  nudge_y = 30, label.r = 0.05, label.size=0.1)+
  scale_fill_manual(values=c(c_cold, c_warm))
#fig1b <- fig1b + annotate(geom = 'text', label = 'Model', x = Inf, y = Inf,
#                          vjust = 1.5, hjust=1.1,fontface ="italic", size=4)

fig1a <- fig1a+ggtitle("Fossil Observation") + xlim(-2,32)+theme(plot.tag = element_text(face = 'bold'))
fig1b <- fig1b+ggtitle("ForamEcoGENIE Model")+xlim(-2,32)+theme(plot.tag = element_text(face = 'bold'))

fig1a<-fig1a+theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))  
fig1b<-fig1b+theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))  
fig1 <- fig1b/fig1a+ plot_annotation(tag_levels = 'a') 
fig1 %>% ggsave(., file="output/fig1.pdf", dpi=400, width=10, height = 6)

###### Species-level thermal niche shift
figs3 <- ggplot()+ 
  geom_point(data=obs_sp_raw, aes(x=SST, y =Abundance, color=age), size=1.5, alpha=0.1) + 
  geom_line(data=obs_sp_smooth,aes(x=model_x, y=model_y, color=age),linewidth=1) +
  facet_wrap(~species,scales = "free_y") + 
  scale_color_manual(values=c(c_cold,c_warm),
                     labels=c("LGM","Pre-industrial")) + 
    labs(x="Sea surface temperature (°C)", y="Abundance (#)")

figs3 <- figs3  + theme(strip.text = element_text(face = "italic"), legend.position = "none")
figs3 <- figs3 + geom_vline(data = thermal_opt(obs_sp_smooth), aes(xintercept = model_x, color = age), linetype = "dashed")

figs3 <- figs3+ geom_label_repel(data = thermal_opt(obs_sp_smooth),
                aes(x = model_x, y=model_y, fill = age, label=round(model_x)),
                color="white",
                nudge_y = 30, label.r = 0.05, label.size=0.1)+
  scale_fill_manual(values=c(c_cold,c_warm))

figs3+ geom_label(data = thermal_opt(obs_sp_smooth),
                        aes(x = model_x, y=model_y, fill = age, label=round(model_x)),
                  position=position_jitter(width=4,height=100))+
  scale_fill_manual(values=c(c_cold,c_warm))

figs3%>%ggsave(file = "output/figs3.jpg", dpi = 400, width = 12, height = 8)

### Fig2
## the same color as the python script
color_palette <- c(c_cold, c_warm, '#420a68', '#932667', '#dd513a', '#fca50a')

genie_fg_smooth$age <- factor(genie_fg_smooth$age, levels = c("lgm","pi","historical","future1p5","future2","future3","future4"))
fig2b <- genie_fg_smooth %>% filter(age!="historical")%>%
  ggplot(aes(x=model_x, y=model_y, color=age))+
  geom_line(linewidth=1) +
  facet_wrap(~species,scales = "free_y", nrow=1) +
  labs(x="Sea surface temperature (°C)",
       y = expression("Abundance (" * "#/m"^3 * ")"),  # LaTeX expression for the y-axis label
       color = "")

fig2b <- fig2b+ 
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8),
        axis.text.y = element_text(size = 8),  # Adjust the font size here (smaller value)
        strip.text = element_text(size = 12),
        strip.background = element_blank()) + 
  scale_color_manual(values = color_palette, 
                     labels = c("Last Glacial Maximum", "Pre-industrial", "2100 (+1°C)", "2100 (+2°C)", "2100 (+3°C)", "2100 (+4°C)"))

ggsave("output/fig2b.jpg", width=9, height=2.5, dpi=300)

## as fig1 but plot chl
figs5 <- ggplot() + geom_line(data=genie_fg_smooth_chl%>% filter(age!="historical"), aes(x=model_x, y =model_y, color=age),  linewidth=1)
figs5 <- figs5 + facet_wrap(~species, scale="free_y",nrow=1)
figs5 <- figs5  + 
  scale_color_manual(values = color_palette, 
                     labels = c("Last Glacial Maximum", "Pre-industrial", "2100 (+1.5°C)", "2100 (+2°C)", "2100 (+3°C)", "2100 (+4°C)"))+
    labs(x=expression("Total Chlorophyll (" * "mg/m"^3 * ")"),
       y = expression("Abundance (" * "#/m"^3 * ")"))
figs5 <- figs5 +  theme(legend.position = "bottom")
ggsave("output/figs5.jpg",figs5, width=8, height=4, dpi=300)
