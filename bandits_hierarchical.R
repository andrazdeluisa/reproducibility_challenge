# Set current folder as working directory
if (rstudioapi::isAvailable()) {
  setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "."))
}

library(ggplot2)

source("utils.R")

regrets_hier <- function(N, L, N_max=3000, seeds){
  # Run HTS algorithm with multiple seeds
  regrets_tsh <- c()

  for (i in 1:length(seeds)){
    tmp <- prepare_tree(N, L, N_max, seeds[i])

    tree <- tmp$tree
    theta <- tmp$theta
    shapes1_t <- tmp$shapes1_t
    shapes2_t <- tmp$shapes2_t
    shapes1_a <- tmp$shapes1_a
    shapes2_a <- tmp$shapes2_a

    regrets_tsh[i] <- tsh(shapes1_t, shapes2_t, shapes1_a, shapes2_a, theta, tree, N_max)
  }
  return(regrets_tsh)
}

set.seed(1)
seeds <- sample(1:1000, 50, replace=FALSE)
N_max <- 3000

# Plot

Ns <- c(50, 100, 1000, 5000)

res_hier <- data.frame(regret=NULL, sd=NULL, N=NULL, L=NULL)

for (k in 1:length(Ns)){
  N <- Ns[k]
  Ls <- c(0, 1, 2, 3, floor(log2(N)))
  for (L in Ls){
    tmp <- regrets_hier(N, L, N_max, seeds)
    res_hier <- rbind(res_hier, data.frame(regret=mean(tmp), sd=sd(tmp), N=as.character(N), L=as.character(L)))
  }
}

res_hier[as.integer(res_hier$L) > 3,]$L <- "log2(N)"

plt <- ggplot(res_hier, aes(x=L)) + geom_point(aes(y=regret, color=N), size=2) +
  geom_errorbar(aes(min=regret-sd, max=regret+sd, color=N), size=1, width=0.25) +
  ylab("Cumulative regret") +
  scale_color_manual(values=c("blue", "orange", "dark green", "red"), breaks=c("50", "100", "1000", "5000"), labels=c("HTS N = 50", "HTS N = 100", "HTS N = 1000", "HTS N = 5000")) +
  ylim(10, 1050) +
  theme(legend.title=element_blank(), legend.position=c(.85, .8), legend.background=element_blank())
ggsave("figures/HTS.pdf", plt, width=5, dpi=600)
