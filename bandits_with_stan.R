# Set current folder as working directory
if (rstudioapi::isAvailable()) {
  setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "."))
}

# Load required libraries
#library(rstan)
library(cmdstanr)
library(posterior)
library(ggplot2)
library(tidyquant)

mod <- cmdstan_model("Models/bernoulli_bandits.stan")

# Simulate bandits
K <- 3
theta <- c(0.6, 0.4, 0.2)

MAX_N <- 200

theta_hat <- matrix(0, MAX_N, K)
p_best <- matrix(0, MAX_N, K)

y <- array(0.0, 0)
z <- array(0.0, 0)
prefix <- function(y, n) array(y, dim = n - 1)
for (n in 1:MAX_N) {
  data <- list(K = K, N = n - 1, y = prefix(y, n), z = prefix(z, n))
  fit <- mod$sample(data, chains = 1, seed=1, iter_warmup=100, iter_sampling=400,
                  init = 0, step_size=0.1, adapt_delta=0.99)
  p_best[n, ] <- fit$summary(c("is_best[1]", "is_best[2]", "is_best[3]"))$mean
  z[n] <- sample(K, 1, replace = TRUE, p_best[n, ])
  y[n] <- rbinom(1, 1, theta[z[n]])
  Sys.sleep(0.2)
}

plt <- ggplot(data.frame(trial = 1:MAX_N, prob_best = p_best[1:MAX_N, 1])) +
  geom_line(aes(trial, prob_best)) +
  geom_ma(aes(trial, prob_best), ma_fun = SMA, n=50)
