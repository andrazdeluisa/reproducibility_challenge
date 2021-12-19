# Set current folder as working directory
if (rstudioapi::isAvailable()) {
  setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "."))
}

library(ggplot2)
library(mvtnorm)

source("utils.R")

compare_contextual_regrets <- function(N, K, eps, d=5, Ts=c(3000), seeds=c(1), random=FALSE){
  # Run contextual bandits algorithms with multiple seeds
  N_max <- max(Ts)
  regrets_lintsc <- matrix(0, nrow=length(seeds), ncol=length(Ts))
  regrets_lints <- matrix(0, nrow=length(seeds), ncol=length(Ts))
  regrets_linucb <- matrix(0, nrow=length(seeds), ncol=length(Ts))
  regrets_linucbc <- matrix(0, nrow=length(seeds), ncol=length(Ts))
  
  for (i in 1:length(seeds)){
    tmp <- prepare_contextual_cluster(N, K, eps, d, seeds[i], random)
    
    clusters <- tmp$clusters
    theta_a <- tmp$theta_a
    
    regrets_lints[i,] <- lin_ts(theta_a, N_max)[Ts]
    regrets_lintsc[i,] <- lin_tsc(clusters, theta_a, N_max)[Ts]
    regrets_linucb[i,] <- lin_ucb(theta_a, N_max)[Ts]
    regrets_linucbc[i,] <- lin_ucbc(clusters, theta_a, N_max)[Ts]
  }
  return(list(lints=regrets_lints, lintsc=regrets_lintsc,
              linucb=regrets_linucb, linucbc=regrets_linucbc))
}

set.seed(1)
seeds <- sample(1:1000, 25, replace=FALSE)
d <- 5
Ts <- 1:20 * 500

# Plot 1

N <- 400
K <- 20
eps <- 0.5

res_1 <- data.frame(regret=0, sd=0, algorithm=c("LinTS", "LinTSC", "LinUCB", "LinUCBC"), t=0)
tmp <- compare_contextual_regrets(N, K, eps, Ts=Ts, seeds=seeds)

for (i in 1:length(Ts)){
  res_1 <- rbind(res_1, data.frame(regret=mean(tmp$lints[,i]), sd=sd(tmp$lints[,i]), algorithm="LinTS", t=Ts[i]))
  res_1 <- rbind(res_1, data.frame(regret=mean(tmp$lintsc[,i]), sd=sd(tmp$lintsc[,i]), algorithm="LinTSC", t=Ts[i]))
  res_1 <- rbind(res_1, data.frame(regret=mean(tmp$linucb[,i]), sd=sd(tmp$linucb[,i]), algorithm="LinUCB", t=Ts[i]))
  res_1 <- rbind(res_1, data.frame(regret=mean(tmp$linucbc[,i]), sd=sd(tmp$linucbc[,i]), algorithm="LinUCBC", t=Ts[i]))
}

plt1 <- ggplot(res_1, aes(x=t)) + geom_point(aes(y=regret, color=algorithm), size=2) +
  geom_line(aes(y=regret, color=algorithm), size=.7) +
  geom_errorbar(aes(min=regret-sd, max=regret+sd, color=algorithm), size=.6, width=200) +
  ylab("Cumulative regret") +
  scale_x_continuous(breaks=c(0, 500 * seq(2, 20, 2)), labels=c(0, 500 * seq(2, 20, 2))) +
  theme(legend.title=element_blank(), legend.position=c(.15, .8), legend.background=element_blank())
ggsave("figures/context1.pdf", plt1, width=5, dpi=600)

# Plot 2

N <- 900
K <- 30
eps <- 0.5

res_2 <- data.frame(regret=0, sd=0, algorithm=c("LinTS", "LinTSC", "LinUCB", "LinUCBC"), t=0)
tmp <- compare_contextual_regrets(N, K, eps, Ts=Ts, seeds=seeds)

for (i in 1:length(Ts)){
  res_2 <- rbind(res_2, data.frame(regret=mean(tmp$lints[,i]), sd=sd(tmp$lints[,i]), algorithm="LinTS", t=Ts[i]))
  res_2 <- rbind(res_2, data.frame(regret=mean(tmp$lintsc[,i]), sd=sd(tmp$lintsc[,i]), algorithm="LinTSC", t=Ts[i]))
  res_2 <- rbind(res_2, data.frame(regret=mean(tmp$linucb[,i]), sd=sd(tmp$linucb[,i]), algorithm="LinUCB", t=Ts[i]))
  res_2 <- rbind(res_2, data.frame(regret=mean(tmp$linucbc[,i]), sd=sd(tmp$linucbc[,i]), algorithm="LinUCBC", t=Ts[i]))
}

plt2 <- ggplot(res_2, aes(x=t)) + geom_point(aes(y=regret, color=algorithm), size=2) +
  geom_line(aes(y=regret, color=algorithm), size=.7) +
  geom_errorbar(aes(min=regret-sd, max=regret+sd, color=algorithm), size=.6, width=200) +
  ylab("Cumulative regret") +
  scale_x_continuous(breaks=c(0, 500 * seq(2, 20, 2)), labels=c(0, 500 * seq(2, 20, 2))) +
  theme(legend.title=element_blank(), legend.position=c(.15, .8), legend.background=element_blank())
ggsave("figures/context2.pdf", plt2, width=5, dpi=600)

# Plot 3

N <- 900
K <- 30
eps <- 0.1

res_3 <- data.frame(regret=0, sd=0, algorithm=c("LinTS", "LinTSC", "LinUCB", "LinUCBC"), t=0)
tmp <- compare_contextual_regrets(N, K, eps, Ts=Ts, seeds=seeds)

for (i in 1:length(Ts)){
  res_3 <- rbind(res_3, data.frame(regret=mean(tmp$lints[,i]), sd=sd(tmp$lints[,i]), algorithm="LinTS", t=Ts[i]))
  res_3 <- rbind(res_3, data.frame(regret=mean(tmp$lintsc[,i]), sd=sd(tmp$lintsc[,i]), algorithm="LinTSC", t=Ts[i]))
  res_3 <- rbind(res_3, data.frame(regret=mean(tmp$linucb[,i]), sd=sd(tmp$linucb[,i]), algorithm="LinUCB", t=Ts[i]))
  res_3 <- rbind(res_3, data.frame(regret=mean(tmp$linucbc[,i]), sd=sd(tmp$linucbc[,i]), algorithm="LinUCBC", t=Ts[i]))
}

plt3 <- ggplot(res_3, aes(x=t)) + geom_point(aes(y=regret, color=algorithm), size=2) +
  geom_line(aes(y=regret, color=algorithm), size=.7) +
  geom_errorbar(aes(min=regret-sd, max=regret+sd, color=algorithm), size=.6, width=200) +
  ylab("Cumulative regret") +
  scale_x_continuous(breaks=c(0, 500 * seq(2, 20, 2)), labels=c(0, 500 * seq(2, 20, 2))) +
  theme(legend.title=element_blank(), legend.position=c(.15, .8), legend.background=element_blank())
ggsave("figures/context3.pdf", plt3, width=5, dpi=600)

### Appendix - Results beyond original paper

# Re-run same experiments with random clustering

N <- 400
K <- 20
eps <- 0.5

res_random <- data.frame(regret=0, sd=0, algorithm=c("LinTS", "LinTSC", "LinUCB", "LinUCBC"), t=0)
tmp <- compare_contextual_regrets(N, K, eps, Ts=Ts, seeds=seeds, random=TRUE)

for (i in 1:length(Ts)){
  res_random <- rbind(res_random, data.frame(regret=mean(tmp$lints[,i]), sd=sd(tmp$lints[,i]), algorithm="LinTS", t=Ts[i]))
  res_random <- rbind(res_random, data.frame(regret=mean(tmp$lintsc[,i]), sd=sd(tmp$lintsc[,i]), algorithm="LinTSC", t=Ts[i]))
  res_random <- rbind(res_random, data.frame(regret=mean(tmp$linucb[,i]), sd=sd(tmp$linucb[,i]), algorithm="LinUCB", t=Ts[i]))
  res_random <- rbind(res_random, data.frame(regret=mean(tmp$linucbc[,i]), sd=sd(tmp$linucbc[,i]), algorithm="LinUCBC", t=Ts[i]))
}

plt_random <- ggplot(res_random, aes(x=t)) + geom_point(aes(y=regret, color=algorithm), size=2) +
  geom_line(aes(y=regret, color=algorithm), size=.7) +
  geom_errorbar(aes(min=regret-sd, max=regret+sd, color=algorithm), size=.6, width=200) +
  ylab("Cumulative regret") +
  scale_x_continuous(breaks=c(0, 500 * seq(2, 20, 2)), labels=c(0, 500 * seq(2, 20, 2))) +
  theme(legend.title=element_blank(), legend.position=c(.15, .8), legend.background=element_blank())
ggsave("figures/context_random.pdf", plt_random, width=5, dpi=600)
