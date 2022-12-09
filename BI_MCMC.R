#1. Packages ####
library(tidyverse)
library(magrittr)
library(lubridate)

#2. Data import ####
load("./DSL_data_2_oclock.RData")

#3. Bayesian inference functions ####
### Notes: use MCMC function on a data.frame (dt) with the following variables (
### OC_c, OC_u, EC_c, EC_u, SIA_c, SIA_u).
MCMC <- function(dt, sps = "S4", k1_avg = 2, k1_sd = 1, k2_avg = 0.4, k2_sd = 0.2){
  dt %<>%
    rename("SIA_c" = paste0(sps, "_c"),
           "SIA_u" = paste0(sps, "_u"))
  
  likelihood <- function(dt, k1, k2){
    dt %>%
      mutate(mean = k1 * EC_c + k2 * SIA_c,
             sd = sqrt(OC_u ^ 2 + (EC_u * k1) ^ 2 + (SIA_u * k2) ^ 2)
      ) %>%
      rename(x = OC_c) %>%
      select(x, mean, sd) %>%
      pmap_dbl(dnorm) %>%
      prod()
  }
  
  Q <- function(dt, k1, k2){
    dnorm(k1, k1_avg, k1_sd) * dnorm(k2, k2_avg, k2_sd) * likelihood(dt, k1, k2)
  }
  
  len <- 5000
  
  k1_0 <- runif(1, 1, 5)
  k2_0 <- runif(1, 0.1, 0.75)
  res <- data.frame(K1 = rep(k1_0, len),
                    K2 = rep(k2_0, len))
  
  for(i in 1:(len - 1)){
    k1_n <- runif(1, 1, 5)
    k2_n <- runif(1, 0.1, 0.75)
    t <- runif(1)
    cond <- (t < Q(dt, k1_n, k2_n)/Q(dt, res$K1[i], res$K2[i]))
    if(cond & (!is.na(cond))){
      res$K1[i+1] <- k1_n
      res$K2[i+1] <- k2_n
    } else {
      res$K1[i+1] <- res$K1[i]
      res$K2[i+1] <- res$K2[i]
    }
    
    if(i%%500 == 0) message(paste0("i = ", i))
  }
  
  message(paste0("Current local time:", Sys.time()))
  
  res %>%
    slice(seq(round(len/5), len, 50)) %>%
    summarise(K1_m = mean(K1),
              K2_m = mean(K2),
              K1_sd = sd(K1),
              K2_sd = sd(K2))
}

MCMC2 <- function(dt, k1_avg = 2, k1_sd = 1, k2_avg = 0.2, k2_sd = 0.1, k3_avg = 0.2, k3_sd = 0.1){
  likelihood <- function(dt, k1, k2, k3){
    dt %>%
      mutate(mean = k1 * EC_c + k2 * S4_c + k3 * N3_c,
             sd = sqrt(OC_u ^ 2 + (EC_u * k1) ^ 2 + (S4_u * k2) ^ 2 + (N3_u * k3) ^ 2)
      ) %>%
      rename(x = OC_c) %>%
      select(x, mean, sd) %>%
      pmap_dbl(dnorm) %>%
      prod()
  }
  
  Q <- function(dt, k1, k2, k3){
    dnorm(k1, k1_avg, k1_sd) * dnorm(k2, k2_avg, k2_sd) * dnorm(k3, k3_avg, k3_sd) * likelihood(dt, k1, k2, k3)
  }
  
  len <- 10000
  
  k1_0 <- runif(1, 1, 5)
  k2_0 <- runif(1, 0.1, 0.55)
  k3_0 <- runif(1, 0.1, 0.55)
  res <- data.frame(K1 = rep(k1_0, len),
                    K2 = rep(k2_0, len),
                    K3 = rep(k3_0, len))
  
  for(i in 1:(len - 1)){
    k1_n <- runif(1, 1, 5)
    k2_n <- runif(1, 0.1, 0.55)
    k3_n <- runif(1, 0.1, 0.55)
    t <- runif(1)
    cond <- (t < Q(dt, k1_n, k2_n, k3_n)/Q(dt, res$K1[i], res$K2[i],  res$K3[i]))
    if(cond & (!is.na(cond))){
      res$K1[i+1] <- k1_n
      res$K2[i+1] <- k2_n
      res$K3[i+1] <- k3_n
    } else {
      res$K1[i+1] <- res$K1[i]
      res$K2[i+1] <- res$K2[i]
      res$K3[i+1] <- res$K3[i]
    }
    
    if(i%%500 == 0) message(paste0("i = ", i))
  }
  
  message(paste0("Current local time:", Sys.time()))
  
  res %>%
    slice(seq(round(len/5), len, 50)) %>%
    summarise(K1_m = mean(K1),
              K2_m = mean(K2),
              K3_m = mean(K3),
              K1_sd = sd(K1),
              K2_sd = sd(K2),
              K3_sd = sd(K3))
}

# 4. Calculation ####
set.seed(1234)
MCMC(DSL_02hr, "S4")


####   END   ####