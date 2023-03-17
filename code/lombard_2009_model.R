library(data.table)
library(tidyverse)

ut1 = c(0.23, 0.3, 0.17, 0.24, 0.1, 0.25, 0.13, 0.37)
t1 = 293 #20 C
ta = c(626, 2155, 7105, 1086, 6876, 6482, 6584, 6584)
# temperature lower bound (in Kel)
tl = c(289.8, 289.3, 284.9, 292.3, 280.1, 281.1, 277.8, NA)
# temperature upper bound (in Kel)
th = c(304, 304.4, 301.8, 303.4, 298.8, 298.5, 293.7, 279.8)
# Arrhenius temperatures
tal = c(30953, 94385, 349991, 53765, 210746, 295374, 300239, NA)
tah = c(249412, 171209,130852,165490,57855,159618,110583,59491)
Species = c("O. universa",
            "G. sacculifer",
            "G. siphonifera",
            "G. ruber w",
            "N. dutertrei",
            "G. bulloides",
            "N. incompta",
            "N. pachyderma"
            )
lombard_data <- data.table(Species, ut1, ta, tl, th,tal,tah) %>% t() %>%
  as.data.table()
names(lombard_data) <- as.character(lombard_data[1,])
lombard_data <- lombard_data[-1, ]
lombard_data <- sapply(lombard_data, function(x) as.numeric(as.character(x))) %>% as.data.frame()
rownames(lombard_data) <- c("ut1", "ta", "tl", "th", "tal", "tah")
lombard_data <- as_tibble(lombard_data)
lombard_data%>%pivot_longer(`O. universa`:`N. pachyderma`, names_to = "Species", values_to = "Parameter_value")

SST = seq(0, 30, 0.1)

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
data.frame(SST, Species="N. dutertrei", Count=growth_rate(SST, lombard_data$`N. dutertrei`)*1000),
data.frame(SST, Species="G. ruber w", Count=growth_rate(SST, lombard_data$`G. ruber w`)*1000),
data.frame(SST, Species="T. sacculifer", Count=growth_rate(SST, lombard_data$`G. sacculifer`)*200),
data.frame(SST, Species="O. universa", Count=growth_rate(SST, lombard_data$`O. universa`)*150),
data.frame(SST, Species="N. dutertrei", Count=growth_rate(SST, lombard_data$`N. dutertrei`)*1000),
data.frame(SST, Species="N. pachyderma", Count=growth_rate(SST, lombard_data$`N. pachyderma`)*1000),
data.frame(SST, Species="N. incompta", Count=growth_rate(SST, lombard_data$`N. incompta`)*4000),
data.frame(SST, Species="G. bulloides", Count=growth_rate(SST, lombard_data$`G. bulloides`)*800),
)
modern_niche$age="modern"
modern_niche$SSS = NA

