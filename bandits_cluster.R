# Set current folder as working directory
if (rstudioapi::isAvailable()) {
  setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "."))
}

library(ggplot2)

source("utils.R")

compare_regrets <- function(N, K, A_, w_, d, N_max=3000, seeds){
  # Run TS algorithms with multiple seeds
  regrets_tsc <- c()
  regrets_ts <- c()
  
  for (i in 1:length(seeds)){
    tmp <- prepare_cluster(N, K, A_, w_, d, N_max, seeds[i])
    
    clusters <- tmp$clusters
    theta <- tmp$theta
    shapes1_c <- tmp$shapes1_c
    shapes2_c <- tmp$shapes2_c
    shapes1_a <- tmp$shapes1_a
    shapes2_a <- tmp$shapes2_a

    regrets_ts[i] <- ts_(shapes1_a, shapes2_a, theta, N_max)
    regrets_tsc[i] <- tsc(shapes1_c, shapes2_c, shapes1_a, shapes2_a, theta, clusters, N_max)
  }
  return(list(ts=regrets_ts, tsc=regrets_tsc))
}

set.seed(1)
seeds <- sample(1:1000, 50, replace=FALSE)
N_max <- 3000

# Plot 1

N <- 100
K <- 10
A_ <- 10
w_ <- .1
ds <- c(.01, .05, .1, .15, .2, .25, .3)

res_d <- data.frame(regret=NULL, sd=NULL, algorithm=NULL)

for (i in 1:length(ds)){
  tmp <- compare_regrets(N, K, A_, w_, ds[i], N_max, seeds)
  res_d <- rbind(res_d, data.frame(regret=mean(tmp$ts), sd=sd(tmp$ts), algorithm="TS", d=ds[i]))
  res_d <- rbind(res_d, data.frame(regret=mean(tmp$tsc), sd=sd(tmp$tsc), algorithm="TSC", d=ds[i]))
}

plt1 <- ggplot(res_d, aes(x=d)) + geom_point(aes(y=regret, color=algorithm), size=3) +
      geom_errorbar(aes(min=regret-sd, max=regret+sd, color=algorithm), size=1, width=0.03) +
      ylab("Cumulative regret") +
      ylim(0, 650) +
      scale_x_continuous(breaks=ds, labels=ds) +
      theme(legend.title=element_blank(), legend.position=c(.85, .5), legend.background=element_blank())
ggsave("figures/TSC_1.pdf", plt1, width=5, dpi=600)

# Plot 2

N <- 100
K <- 10
A_ <- 10
w_ <- c(0, .1, .2, .3)
d <- .1

res_w <- data.frame(regret=NULL, sd=NULL, algorithm=NULL)

for (i in 1:length(w_)){
  tmp <- compare_regrets(N, K, A_, w_[i], d, N_max, seeds)
  res_w <- rbind(res_w, data.frame(regret=mean(tmp$ts), sd=sd(tmp$ts), algorithm="TS", w=w_[i]))
  res_w <- rbind(res_w, data.frame(regret=mean(tmp$tsc), sd=sd(tmp$tsc), algorithm="TSC", w=w_[i]))
}

plt2 <- ggplot(res_w, aes(x=w)) + geom_point(aes(y=regret, color=algorithm), size=3) +
  geom_errorbar(aes(min=regret-sd, max=regret+sd, color=algorithm), size=1, width=.03) +
  ylab("Cumulative regret") + xlab("w*") +
  ylim(0, 650) +
  scale_x_continuous(breaks=w_, labels=w_) +
  theme(legend.title=element_blank(), legend.position=c(.85, .5), legend.background=element_blank())
ggsave("figures/TSC_2.pdf", plt2, width=5, dpi=600)

# Plot 3

N <- c(10, 30, 50, 70, 90, 110, 130, 150)
K <- floor(sqrt(N))
A_ <- floor(sqrt(N))
w_ <- .1
d <- .1

res_n <- data.frame(regret=NULL, sd=NULL, algorithm=NULL)

for (i in 1:length(N)){
  tmp <- compare_regrets(N[i], K[i], A_[i], w_, d, N_max, seeds)
  res_n <- rbind(res_n, data.frame(regret=mean(tmp$ts), sd=sd(tmp$ts), algorithm="TS", n=N[i]))
  res_n <- rbind(res_n, data.frame(regret=mean(tmp$tsc), sd=sd(tmp$tsc), algorithm="TSC", n=N[i]))
}

plt3 <- ggplot(res_n, aes(x=n)) + geom_point(aes(y=regret, color=algorithm), size=3) +
  geom_errorbar(aes(min=regret-sd, max=regret+sd, color=algorithm), size=1, width=10) +
  ylab("Cumulative regret") + xlab("N") +
  ylim(0, 650) +
  scale_x_continuous(breaks=N, labels=N) +
  theme(legend.title=element_blank(), legend.position=c(.85, .55), legend.background=element_blank())
ggsave("figures/TSC_3.pdf", plt3, width=5, dpi=600)

# Plot 4

N <- 100
K <- c(2, 5, 10, 15, 20, 25)
A_ <- 10
w_ <- .1
d <- .1

res_k <- data.frame(regret=NULL, sd=NULL, algorithm=NULL)

for (i in 1:length(K)){
  tmp <- compare_regrets(N, K[i], A_, w_, d, N_max, seeds)
  res_k <- rbind(res_k, data.frame(regret=mean(tmp$ts), sd=sd(tmp$ts), algorithm="TS", k=K[i]))
  res_k <- rbind(res_k, data.frame(regret=mean(tmp$tsc), sd=sd(tmp$tsc), algorithm="TSC", k=K[i]))
}

plt4 <- ggplot(res_k, aes(x=k)) + geom_point(aes(y=regret, color=algorithm), size=3) +
  geom_errorbar(aes(min=regret-sd, max=regret+sd, color=algorithm), size=1, width=2) +
  ylab("Cumulative regret") + xlab("K") +
  ylim(0, 650) +
  scale_x_continuous(breaks=K, labels=K) +
  theme(legend.title=element_blank(), legend.position=c(.85, .2), legend.background=element_blank())
ggsave("figures/TSC_4.pdf", plt4, width=5, dpi=600)

# Plot 5

N <- 100
K <- 10
A_ <- 10 * 1:5
w_ <- .1
d <- .1

res_a <- data.frame(regret=NULL, sd=NULL, algorithm=NULL)

for (i in 1:length(A_)){
  tmp <- compare_regrets(N, K, A_[i], w_, d, N_max, seeds)
  res_a <- rbind(res_a, data.frame(regret=mean(tmp$ts), sd=sd(tmp$ts), algorithm="TS", a=A_[i]))
  res_a <- rbind(res_a, data.frame(regret=mean(tmp$tsc), sd=sd(tmp$tsc), algorithm="TSC", a=A_[i]))
}

plt5 <- ggplot(res_a, aes(x=a)) + geom_point(aes(y=regret, color=algorithm), size=3) +
  geom_errorbar(aes(min=regret-sd, max=regret+sd, color=algorithm), size=1, width=4) +
  ylab("Cumulative regret") + xlab("A*") +
  ylim(0, 650) +
  scale_x_continuous(breaks=A_, labels=A_) +
  theme(legend.title=element_blank(), legend.position=c(.85, .8), legend.background=element_blank())
ggsave("figures/TSC_5.pdf", plt5, width=5, dpi=600)
