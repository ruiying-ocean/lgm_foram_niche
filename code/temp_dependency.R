library(tidyverse)

sst <- seq(0,30,1)
eppley <- 0.59 * exp(0.0633*sst)
medusa <- 1.066**sst
bissinger <-  0.81 * exp(0.0631*sst)
temp_a <- 0.05
temp_T0 <- 20
ward <- exp(temp_a*(sst-temp_T0))
d <- data.frame(sst = sst, eppley = eppley,
           medusa=medusa, bissinger=bissinger, ward=ward)
d <- d |> tidyr::pivot_longer(cols=eppley:ward, names_to = "scheme", values_to = "gamma")
ggplot(d, aes(x=sst, y=gamma)) + geom_line(aes(color=scheme))

# Q10 value = (gamma2/gamma1)^（10/(T2-T1))
d %>% group_by(scheme) %>% summarise(q10 = (max(gamma)/min(gamma))**1/3)

# Data: Andy Fraass 2015 annual review
extinct.rate  <- c(NA, 0.324, 0.445, 0.265, 0.499, 0.510 , 0.407,0.971, NA, 0.353, NA, 17.028, 0.968, 0.433, 0.271)
origin.rate <- c(0.313, NA, NA, 0.458, 0.849, 0.510, 0.538, 1.527, 0.374, NA, 10.674, 36.620, NA, 0.411, NA)
ggplot(data.frame(E=extinct.rate, O=origin.rate), aes(x=E, y=O))+ geom_point() + geom_smooth(method="lm")
