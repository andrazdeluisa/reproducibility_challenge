# Set current folder as working directory
if (rstudioapi::isAvailable()) {
  setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "."))
}

library(ggplot2)

source("utils.R")

compare_regrets_violated_hier <- function(N, K, L, Ts, seeds){
  # Run hierarchical TS algorithms with multiple seeds
  
  N_max <- max(Ts)
  regrets_hts <- matrix(0, nrow=length(seeds), ncol=length(Ts))
  regrets_uct <- matrix(0, nrow=length(seeds), ncol=length(Ts))
  
  for (i in 1:length(seeds)){
    tmp <- prepare_tree_violated(N, K, L, seeds[i])
    
    tree <- tmp$tree
    theta <- tmp$theta
    
    shapes1_a <- tmp$shapes1_a
    shapes2_a <- tmp$shapes2_a
    shapes1_t <- tmp$shapes1_t
    shapes2_t <- tmp$shapes2_t
    
    regrets_hts[i,] <- hts_non_binary(shapes1_t, shapes2_t, shapes1_a, shapes2_a, theta, tree, N_max)[Ts]
    regrets_uct[i,] <- uct(shapes1_t, shapes1_a, theta, tree, N_max)[Ts]
  }
  return(list(hts=regrets_hts, uct=regrets_uct))
}


### Plot

N <- 5000
K <- 15
set.seed(1)
seeds <- sample(1:1000, size=100, replace=FALSE)
Ts <- 1:10 * 500

res <- data.frame(regret=0, sd=0, algorithm=c("TSC", "HTS L = 2", "HTS L = 3", "UCT L = 2", "UCT L = 3"), t=0)

tmp <- compare_regrets_violated_hier(N, K, L=1, Ts=Ts, seeds=seeds)
for (i in 1:length(Ts)){
  res <- rbind(res, data.frame(regret=mean(tmp$hts[,i]), sd=sd(tmp$hts[,i]), algorithm="TSC", t=Ts[i]))
}

tmp <- compare_regrets_violated_hier(N, K, L=2, Ts=Ts, seeds=seeds)
for (i in 1:length(Ts)){
  res <- rbind(res, data.frame(regret=mean(tmp$hts[,i]), sd=sd(tmp$hts[,i]), algorithm="HTS L = 2", t=Ts[i]))
  res <- rbind(res, data.frame(regret=mean(tmp$uct[,i]), sd=sd(tmp$uct[,i]), algorithm="UCT L = 2", t=Ts[i]))
}

tmp <- compare_regrets_violated_hier(N, K, L=3, Ts=Ts, seeds=seeds)
for (i in 1:length(Ts)){
  res <- rbind(res, data.frame(regret=mean(tmp$hts[,i]), sd=sd(tmp$hts[,i]), algorithm="HTS L = 3", t=Ts[i]))
  res <- rbind(res, data.frame(regret=mean(tmp$uct[,i]), sd=sd(tmp$uct[,i]), algorithm="UCT L = 3", t=Ts[i]))
}

plt <- ggplot(res, aes(x=t)) + geom_point(aes(y=regret, color=algorithm), size=2) +
  geom_line(aes(y=regret, color=algorithm), size=.7) +
  geom_errorbar(aes(min=regret-sd, max=regret+sd, color=algorithm), size=.6, width=200) +
  ylab("Cumulative regret") +
  scale_x_continuous(breaks=c(0, 500 * seq(2, 10, 2)), labels=c(0, 500 * seq(2, 10, 2))) +
  theme(legend.title=element_blank(), legend.position=c(.15, .8), legend.background=element_blank())
ggsave("figures/violated_hier.pdf", plt, width=5, height=5, dpi=600)

### Plot 2 - results beyond original paper

N <- 1000
K <- 10
set.seed(1)
seeds <- sample(1:1000, size=100, replace=FALSE)
Ts <- 1:10 * 1000

res_2 <- data.frame(regret=0, sd=0, algorithm=c("TSC", "HTS L = 2", "UCT L = 2"), t=0)

tmp <- compare_regrets_violated_hier(N, K, L=1, Ts=Ts, seeds=seeds)
for (i in 1:length(Ts)){
  res_2 <- rbind(res_2, data.frame(regret=mean(tmp$hts[,i]), sd=sd(tmp$hts[,i]), algorithm="TSC", t=Ts[i]))
}

tmp <- compare_regrets_violated_hier(N, K, L=2, Ts=Ts, seeds=seeds)
for (i in 1:length(Ts)){
  res_2 <- rbind(res_2, data.frame(regret=mean(tmp$hts[,i]), sd=sd(tmp$hts[,i]), algorithm="HTS L = 2", t=Ts[i]))
  res_2 <- rbind(res_2, data.frame(regret=mean(tmp$uct[,i]), sd=sd(tmp$uct[,i]), algorithm="UCT L = 2", t=Ts[i]))
}

plt_2 <- ggplot(res_2, aes(x=t)) + geom_point(aes(y=regret, color=algorithm), size=2) +
  geom_line(aes(y=regret, color=algorithm), size=.7) +
  geom_errorbar(aes(min=regret-sd, max=regret+sd, color=algorithm), size=.6, width=200) +
  ylab("Cumulative regret") +
  scale_x_continuous(breaks=c(0, 1000 * seq(2, 10, 2)), labels=c(0, 1000 * seq(2, 10, 2))) +
  theme(legend.title=element_blank(), legend.position=c(.15, .8), legend.background=element_blank())
ggsave("figures/violated_hier_2.pdf", plt_2, width=5, height=5, dpi=600)
