library(readr)
library(ggfortify)

## using delta abundance for multivariate linear model

df <- read_csv("data/delta_4pca.csv")
df_bn <- df[,-c(1,3,4,5)]
df_bs <- df[,-c(1,2,4,5)]
df_sn <- df[,-c(1,2,3,5)]
df_ss <- df[,-c(1,2,3,4)]
names(df_bn)[1] <- "dabundance"
names(df_bs)[1] <- "dabundance"
names(df_sn)[1] <- "dabundance"
names(df_ss)[1] <- "dabundance"

df_bn <- as_tibble(scale(df_bn))
df_bs <- as_tibble(scale(df_sn))
df_sn <- as_tibble(scale(df_sn))
df_ss <- as_tibble(scale(df_ss))

                                        # ladder transformation -> linear model
lambda <- 0.7
df_bn$new <-  (df_bn$dabundance^ lambda - 1) / lambda
mod_bn  <- lm(new ~ dsst + dsal + dchl + dpo4 + dfe,
              data = df_bn)
autoplot(mod_bn)
summary(mod_bn)

                                        #assess.lm <- gvlma(mod_bn)
                                        #pass the global stat
df_bs$new <-  (df_bs$dabundance^ lambda - 1) / lambda
mod_bs  <- lm(new ~ dsst + dsal + dchl + dpo4 + dfe,
              data = df_bs)
autoplot(mod_bs)
summary(mod_bs)

df_sn$new <-  (df_sn$dabundance^ lambda - 1) / lambda
mod_sn  <- lm(new ~ dsst + dsal + dchl + dpo4 + dfe,
              data = df_sn)
autoplot(mod_sn)
summary(mod_sn)

df_ss$new <-  (df_ss$dabundance^ lambda - 1) / lambda
mod_ss  <- lm(new ~ dsst + dsal + dchl + dpo4 + dfe,
              data = df_ss)
autoplot(mod_ss)
summary(mod_ss)
