# Set current folder as working directory
if (rstudioapi::isAvailable()) {
  setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "."))
}

library(ggplot2)

source("utils.R")

compare_regrets_violated <- function(N, K, N_max=3000, seeds){
  # Run TS algorithms with multiple seeds
  regrets_tsc <- c()
  regrets_ts <- c()
  regrets_ucb1 <- c()
  regrets_tsmax <- c()
  regrets_ucbc <- c()
  
  for (i in 1:length(seeds)){
    tmp <- prepare_kmean_cluster(N, K, seeds[i])
    
    clusters <- tmp$clusters
    theta <- tmp$theta
    shapes1_c <- tmp$shapes1_c
    shapes2_c <- tmp$shapes2_c
    shapes1_a <- tmp$shapes1_a
    shapes2_a <- tmp$shapes2_a

    regrets_ts[i] <- ts_(shapes1_a, shapes2_a, theta, N_max)
    regrets_tsc[i] <- tsc(shapes1_c, shapes2_c, shapes1_a, shapes2_a, theta, clusters, N_max)
    regrets_ucb1[i] <- ucb1(theta, N_max)
    regrets_tsmax[i] <- ts_max(shapes1_c, shapes2_c, shapes1_a, shapes2_a, theta, clusters, N_max)
    regrets_ucbc[i] <- ucbc(shapes1_a, shapes2_a, theta, clusters, N_max)
  }
  return(list(ts=regrets_ts, tsc=regrets_tsc, ucb1=regrets_ucb1, tsmax=regrets_tsmax, ucbc=regrets_ucbc))
}

# Plot 1

N <- 100
K <- 10
Ns <- 1:12 * 250
set.seed(1)
seeds <- sample(1:1000, 100, replace=FALSE)

res_small <- data.frame(regret=0, sd=0, algorithm=c("TS", "TSC", "UCB", "TS Max", "UCBC"), t=0)

for (i in 1:length(Ns)){
  tmp <- compare_regrets_violated(N, K, Ns[i], seeds)
  res_small <- rbind(res_small, data.frame(regret=mean(tmp$ts), sd=sd(tmp$ts), algorithm="TS", t=Ns[i]))
  res_small <- rbind(res_small, data.frame(regret=mean(tmp$tsc), sd=sd(tmp$tsc), algorithm="TSC", t=Ns[i]))
  res_small <- rbind(res_small, data.frame(regret=mean(tmp$ucb1), sd=sd(tmp$ucb1), algorithm="UCB", t=Ns[i]))
  res_small <- rbind(res_small, data.frame(regret=mean(tmp$tsmax), sd=sd(tmp$tsmax), algorithm="TS Max", t=Ns[i]))
  res_small <- rbind(res_small, data.frame(regret=mean(tmp$ucbc), sd=sd(tmp$ucbc), algorithm="UCBC", t=Ns[i]))
}

plt1 <- ggplot(res_small, aes(x=t)) + geom_point(aes(y=regret, color=algorithm), size=2) +
  geom_line(aes(y=regret, color=algorithm), size=.7) +
  geom_errorbar(aes(min=regret-sd, max=regret+sd, color=algorithm), size=.5, width=60) +
  ylab("Cumulative regret") +
  ylim(0, 1000) +
  scale_x_continuous(breaks=c(0, Ns), labels=c(0, Ns)) +
  theme(legend.title=element_blank(), legend.position=c(.15, .75), legend.background=element_blank())
ggsave("figures/v_ass1.pdf", plt1, width=5, dpi=600)

# Plot 2

N <- 1000
K <- 32
Ns <- 1:12 * 250
set.seed(1)
seeds <- sample(1:1000, 100, replace=FALSE)

res_large <- data.frame(regret=0, sd=0, algorithm=c("TS", "TSC", "UCB", "TS Max", "UCBC"), t=0)

for (i in 1:length(Ns)){
  tmp <- compare_regrets_violated(N, K, Ns[i], seeds)
  res_large <- rbind(res_large, data.frame(regret=mean(tmp$ts), sd=sd(tmp$ts), algorithm="TS", t=Ns[i]))
  res_large <- rbind(res_large, data.frame(regret=mean(tmp$tsc), sd=sd(tmp$tsc), algorithm="TSC", t=Ns[i]))
  res_large <- rbind(res_large, data.frame(regret=mean(tmp$ucb1), sd=sd(tmp$ucb1), algorithm="UCB", t=Ns[i]))
  res_large <- rbind(res_large, data.frame(regret=mean(tmp$tsmax), sd=sd(tmp$tsmax), algorithm="TS Max", t=Ns[i]))
  res_large <- rbind(res_large, data.frame(regret=mean(tmp$ucbc), sd=sd(tmp$ucbc), algorithm="UCBC", t=Ns[i]))
}

plt2 <- ggplot(res_large, aes(x=t)) + geom_point(aes(y=regret, color=algorithm), size=2) +
  geom_line(aes(y=regret, color=algorithm), size=.7) +
  geom_errorbar(aes(min=regret-sd, max=regret+sd, color=algorithm), size=.5, width=60) +
  ylab("Cumulative regret") +
  ylim(0, 1290) +
  scale_x_continuous(breaks=c(0, Ns), labels=c(0, Ns)) +
  theme(legend.title=element_blank(), legend.position=c(.15, .75), legend.background=element_blank())
ggsave("figures/v_ass2.pdf", plt2, width=5, dpi=600)
