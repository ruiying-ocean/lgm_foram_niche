library(data.table)
library(tidyverse)

Species = c("O. universa",
            "G. sacculifer",
            "G. siphonifera",
            "G. ruber albus",
            "N. dutertrei",
            "G. bulloides",
            "N. incompta",
            "N. pachyderma"
)

ut1 = c(0.23, 0.30, 0.17, 0.24, 0.10, 0.25, 0.13, 0.37)
ta = c(626, 2155, 7105, 1086, 6876, 6482, 6584, 6584)
# temperature lower bound (in Kel)
tl = c(289.8, 289.3, 284.9, 292.3, 280.1, 281.1, 277.8, 277.8 )
# temperature upper bound (in Kel)
th = c(304, 304.4, 301.8, 303.4, 298.8, 298.5, 293.7, 279.8)
# Arrhenius temperatures
tal = c(30953, 94385, 349991, 53765, 210746, 295374, 300239, 300239)
tah = c(249412, 171209,130852,165490,57855,159618,110583,59491)
lombard_lab <- data.table(Species, ut1, ta, tl, th,tal,tah) %>% t() %>%
  as.data.table()


ut1 = c(0.21, 0.29, 0.19, 0.19, 0.17, 0.31, 0.18, 0.19)
kn = c(1.73, 1.32, 1.19, 0.51, 1.00, 6.84, 3.33, 4.70)
tah = c(74, 313, 102000, 39284, 32319, 52575, 51836, 23802)
th = c(305, 305, 302, 303, 304, 299, 296, 281)
tal = c(31002, 51870, 270000, 44807, 103000, 202000, 164000, 20900)
tl = c(287, 289, 285, 291, 281, 281, 277, 260)
ta = c(5598, 3523, 10427, 7852, 8536, 9006, 8347, 3287)
lombard_foramclim <- data.table(Species, ut1, ta, tl, th,tal,tah,kn) %>% t() %>%
  as.data.table()

names(lombard_lab) <- as.character(lombard_lab[1,])
lombard_lab <- lombard_lab[-1, ]
lombard_lab <- sapply(lombard_lab, function(x) as.numeric(as.character(x))) %>% as.data.frame()
rownames(lombard_lab) <- c("ut1", "ta", "tl", "th", "tal", "tah")
lombard_lab <- as_tibble(lombard_lab)
lombard_lab%>%pivot_longer(`O. universa`:`N. pachyderma`, names_to = "Species", values_to = "Parameter_value")

SST = seq(0, 35, 0.1)

growth_rate <- function(temperature, pars, t1=293) {
  temperature = temperature+273.5
  ut1 = pars[1]
  ta = pars[2]
  tl=pars[3]
  th=pars[4]
  tal= pars[5]
  tah=pars[6]
  fra1 = ut1 * exp(ta/t1 - ta/temperature)
  fra2 = 1+exp(tal/temperature - tal/tl)+exp(tah/th-tah/temperature)
  ut= fra1/fra2
  return(ut)
}

modern_niche <- bind_rows(
data.frame(SST, Species="N. dutertrei", u=growth_rate(SST, lombard_foramclim$`N. dutertrei`)),
data.frame(SST, Species="G. ruber albus", u=growth_rate(SST, lombard_foramclim$`G. ruber albus`)),
data.frame(SST, Species="T. sacculifer", u=growth_rate(SST, lombard_foramclim$`G. sacculifer`)),
data.frame(SST, Species="O. universa", u=growth_rate(SST, lombard_foramclim$`O. universa`)),
data.frame(SST, Species="N. dutertrei", u=growth_rate(SST, lombard_foramclim$`N. dutertrei`)),
data.frame(SST, Species="N. pachyderma", u=growth_rate(SST, lombard_foramclim$`N. pachyderma`)),
data.frame(SST, Species="N. incompta", u=growth_rate(SST, lombard_foramclim$`N. incompta`)),
data.frame(SST, Species="G. bulloides", u=growth_rate(SST, lombard_foramclim$`G. bulloides`)),
)

modern_niche <- modern_niche %>%  mutate(u = if_else(Species == "N. pachyderma" & SST < 5, NA, u))
modern_niche$age="modern"
modern_niche %>% ggplot(aes(x=SST,y=u)) + geom_line(aes(color=Species), linewidth=1) + ylim(c(0,0.5))
