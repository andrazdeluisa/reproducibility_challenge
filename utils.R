# Set current folder as working directory
if (rstudioapi::isAvailable()) {
  setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "."))
}

### Simple clustering

prepare_cluster <- function(N=100, K=10, A_=10, w_=.1, d=.1, N_max=3000, seed=1){
  # Assign clusters and probabilities to each of the arms
  # Prepare data to be processed with Thompson Sampling algorithms
  set.seed(seed)
  
  if (K > 2){
    clusters <- sample(2:K, N, replace=TRUE) # uniformly randomly assign arms to clusters
  } else {
    clusters <- rep(K, N)
  }
  clusters[1:A_] <- 1
  clusters[(A_+1):(A_+K-1)] <- 2:K
  
  theta <- rep(0, N)
  theta[1] <- 0.6
  theta[2] <- 0.6 - w_
  if (A_ > 2){
    theta[3:A_] <- runif(A_-2, min=0.6-w_, max=0.6)
  }
  theta[(A_+1):N] <- runif(N-A_, min=0.5-w_-d, max=0.6-w_-d) # assign theta to each arm, uniformly randomly
  
  for (i in 2:K){
    idx <- which(clusters == i)
    theta[idx[1]] <- 0.6 - w_ - d # max value in each cluster
    if (length(idx) >= 2){
      theta[idx[2]] <- 0.5 - w_ - d # min value in each cluster
    }
  }
  
  # Uniform prior beliefs
  shapes1_c <- rep(1, K)
  shapes2_c <- rep(1, K)
  
  shapes1_a <- rep(1, N)
  shapes2_a <- rep(1, N)
  
  return(list(clusters=clusters, theta=theta,
              shapes1_c=shapes1_c, shapes2_c=shapes2_c, 
              shapes1_a=shapes1_a, shapes2_a=shapes2_a))
}

ts_ <- function(shapes1, shapes2, theta, N_max=3000){
  # Conventional Thompson sampling
  
  rewards <- c()
  N <- length(shapes1)
  for (i in 1:N_max){
    samp <- rbeta(N, shapes1, shapes2)
    exp <- shapes1 / (shapes1 + shapes2)
    k <- order(samp, decreasing=TRUE)[1]
    r <- rbinom(1, 1, prob=theta[k])
    shapes1[k] <- shapes1[k] + r
    shapes2[k] <- shapes2[k] + (1 - r)
    rewards[i] <- r
  }
  regret <- max(theta) * N_max - sum(rewards)
  return(regret)
}

tsc <- function(shapes1_c, shapes2_c, shapes1_a, shapes2_a, theta, clusters, N_max=3000){
  # Thompson sampling algorithm with clustered arms
  
  rewards <- c()
  K <- length(shapes1_c)
  for (i in 1:N_max){
    samp_c <- rbeta(K, shapes1_c, shapes2_c) # sample from cluster posteriors
    k <- order(samp_c, decreasing=TRUE)[1] # select best cluster
    idx <- which(clusters == k) # assume non-empty (by cluster construction)
    samp_a <- rbeta(length(idx), shapes1_a[idx], shapes2_a[idx])
    a <- order(samp_a, decreasing=TRUE)[1] # select best arm in cluster
    r <- rbinom(1, 1, prob=theta[idx[a]])
    shapes1_c[k] <- shapes1_c[k] + r
    shapes2_c[k] <- shapes2_c[k] + (1 - r)
    shapes1_a[idx[a]] <- shapes1_a[idx[a]] + r
    shapes2_a[idx[a]] <- shapes2_a[idx[a]] + (1 - r)
    rewards[i] <- r
  }
  regret <- max(theta) * N_max - sum(rewards)
  return(regret)
}

### Simple clustering with violated assumptions

prepare_kmean_cluster <- function(N=100, K=10, seed=1){
  # Cluster uniformly randomly sampled arms with kmean algorithm
  # Apply a smooth and periodic function to the arms
  # Prepare data with violated assumptions to be processed with Thompson Sampling algorithms
  set.seed(seed)
  
  f <- function(x){(sin(13 * x) * sin(27 * x) + 1) / 2}
  xs <- runif(N, 0, 1)
  theta <- f(xs)
  clusters <- kmeans(xs, K)$cluster
  
  # Uniform prior beliefs
  shapes1_c <- rep(1, K)
  shapes2_c <- rep(1, K)
  
  shapes1_a <- rep(1, N)
  shapes2_a <- rep(1, N)
  
  return(list(clusters=clusters, theta=theta,
              shapes1_c=shapes1_c, shapes2_c=shapes2_c, 
              shapes1_a=shapes1_a, shapes2_a=shapes2_a))
}

ucb1 <- function(theta, N_max=3000){
  # UCB1 deterministic policy (Auer 2002)
  
  rewards <- c()
  N <- length(theta)
  nr_plays <- rep(1, N)
  arm_r <- rbinom(N, 1, theta) # initialization - play each arm once
  rewards[1:N] <- arm_r
  for (i in (N + 1):N_max){
    tmp <- arm_r / nr_plays + sqrt(2 * log(i) / nr_plays) # play arm that maximizes this expression
    a <- order(tmp, decreasing=TRUE)[1]
    r <- rbinom(1, 1, theta[a])
    rewards[i] <- r
    arm_r[a] <- arm_r[a] + r # arm a reward obtained so far
    nr_plays[a] <- nr_plays[a] + 1 # nr of times arm a was played so far
  }
  
  regret <- max(theta) * N_max - sum(rewards[1:N_max])
  return(regret)
}

ts_max <- function(shapes1_c, shapes2_c, shapes1_a, shapes2_a, theta, clusters, N_max=3000){
  # TS Max algorithm (Zhao 2019) for bandits with clustered arms
  
  K <- length(shapes1_c)
  rewards <- c()
  for (i in 1:N_max){
    samp_c <- rbeta(K, shapes1_c, shapes2_c) # sample from cluster posteriors
    k <- order(samp_c, decreasing=TRUE)[1] # select best cluster
    idx <- which(clusters == k) # assume non-empty (by cluster construction)
    samp_a <- rbeta(length(idx), shapes1_a[idx], shapes2_a[idx])
    a <- order(samp_a, decreasing=TRUE)[1] # select best arm in cluster
    r <- rbinom(1, 1, prob=theta[idx[a]])
    shapes1_a[idx[a]] <- shapes1_a[idx[a]] + r
    shapes2_a[idx[a]] <- shapes2_a[idx[a]] + (1 - r)
    arm_max <- idx[order(shapes1_a[idx] / (shapes1_a[idx] + shapes2_a[idx]), decreasing=TRUE)[1]] # maximum expected value
    shapes1_c[k] <- shapes1_a[arm_max] + 1
    shapes2_c[k] <- shapes2_a[arm_max] + 1
    rewards[i] <- r
  }
  regret <- max(theta) * N_max - sum(rewards)
  return(regret)
}

ucbc <- function(shapes1_a, shapes2_a, theta, clusters, N_max=3000){
  # UCBC algorithm (called TLP in Pandey 2007) for bandits with clustered arms
  # With MAX strategy (other options include MEAN and PMAX) based on UCB1 (special case of UCT)
  K <- length(unique(clusters))
  rewards <- c()
  max_arms <- c()
  for (i in 1:K){max_arms[i] <- which(clusters == i)[1]}
  nr_plays_c <- rep(1, K)
  nr_plays_a <- rep(1, N)
  
  for (i in 1:N_max){
    alpha <- shapes1_a[max_arms]
    beta <- shapes2_a[max_arms]
    exp_c <- alpha / (alpha + beta) # max expected theta values in each cluster
    tmp <- exp_c + sqrt(2 * log(i) / nr_plays_c) # play cluster that maximizes this expression
    k <- order(tmp, decreasing=TRUE)[1] # select best cluster to play
    idx <- which(clusters == k) # assume non-empty (by cluster construction)
    
    alpha_a <- shapes1_a[idx]
    beta_a <- shapes2_a[idx]
    exp_a <- alpha_a / (alpha_a + beta_a)
    tmp_a <- exp_a + sqrt(2 * log(nr_plays_c[k]) / nr_plays_a[idx]) # play arm that maximizes this expression
    a <- order(tmp_a, decreasing=TRUE)[1] # select best arm in cluster
    nr_plays_a[idx[a]] <- nr_plays_a[idx[a]] + 1
    nr_plays_c[k] <- nr_plays_c[k] + 1
    
    r <- rbinom(1, 1, prob=theta[idx[a]]) # draw the reward
    shapes1_a[idx[a]] <- shapes1_a[idx[a]] + r
    shapes2_a[idx[a]] <- shapes2_a[idx[a]] + (1 - r)
    
    max_arms[k] <- idx[order(shapes1_a[idx] / (shapes1_a[idx] + shapes2_a[idx]), decreasing=TRUE)[1]] # update arm which represents the cluster
    rewards[i] <- r
  }
  
  regret <- max(theta) * N_max - sum(rewards)
  return(regret)
}

### Hierarchical clustering

prepare_tree <- function(N=100, L=0, N_max=3000, seed=1){
  # Assign probabilities to each of the arms
  # Build a strongly dominant hierarchical tree
  # Prepare data to be processed with Thompson Sampling algorithms
  set.seed(seed)
  
  theta <- runif(N, min=0.1, max=0.8)
  ord <- order(theta, decreasing=TRUE)
  tree <- list(list(ord))
  shapes1_t <- list(1)
  shapes2_t <- list(1)
  shapes1_a <- rep(1, N)
  shapes2_a <- rep(1, N)
  
  if (L > 0) {
    for (i in 2:(L+1)){
      tree[[i]] <- list()
      shapes1_t[[i]] <- rep(1, 2 * 2 ** (i - 2))
      shapes2_t[[i]] <- rep(1, 2 * 2 ** (i - 2))
      for (j in 1:(2 ** (i - 2))){
        node <- tree[[i-1]][[j]]
        len <- length(node)
        left <- node[1:ceiling(len / 2)]
        right <- node[(ceiling(len / 2) + 1):len]
        tree[[i]][[2 * j - 1]] <- left
        tree[[i]][[2 * j]] <- right
      }
    }
  }
  
  return(list(tree=tree, theta=theta,
              shapes1_t=shapes1_t, shapes2_t=shapes2_t, 
              shapes1_a=shapes1_a, shapes2_a=shapes2_a))
}

tsh <- function(shapes1_t, shapes2_t, shapes1_a, shapes2_a, theta, tree, N_max=3000){
  # Thompson sampling with hierarchically clustered arms
  
  rewards <- c()
  L <- length(tree) - 1
  for (i in 1:N_max){
    path <- c(1)
    if (L > 0) {
      for (j in 2:(L + 1)){
        idxs <- (2 * path[length(path)] - 1):(2*path[length(path)])
        node <- tree[[j]][idxs]
        samp <- rbeta(length(node), shapes1_t[[j]][idxs], shapes2_t[[j]][idxs])
        path <- append(path, idxs[order(samp, decreasing=TRUE)[1]])
      }
    }
    node <- tree[[length(tree)]][[path[length(path)]]]
    samp <- rbeta(length(node), shapes1_a[node], shapes2_a[node])
    sel_arm <- node[order(samp, decreasing=TRUE)[1]]
    r <- rbinom(1, 1, theta[sel_arm])
    shapes1_a[sel_arm] <- shapes1_a[sel_arm] + r
    shapes2_a[sel_arm] <- shapes2_a[sel_arm] + (1 - r)
    for (j in 1:(L + 1)){
      shapes1_t[[j]][path[j]] <- shapes1_t[[j]][path[j]] + r
      shapes2_t[[j]][path[j]] <- shapes2_t[[j]][path[j]] + (1 - r)
    }
    rewards[i] <- r
  }
  
  regret <- N_max * max(theta) - sum(rewards)
  return(regret)
}

### Hierarchical clustering with violated assumptions

prepare_tree_violated <- function(N, K, L, seed){
  # Cluster uniformly randomly sampled arms with kmean algorithm recursively
  # Build a hierarchical tree
  # Apply a smooth and periodic function to the arms
  # Prepare data with violated assumptions to be processed with Hierarchical Thompson Sampling algorithms
  set.seed(seed)
  
  f <- function(x){(sin(13 * x) * sin(27 * x) + 1) / 2}
  xs <- runif(N, 0, 1)
  theta <- f(xs)
  cluster_level1 <- kmeans(xs, K)$cluster
  cluster_level2 <- list()
  cluster_level3 <- list()
  
  tree <- list()
  shapes1_t <- list(rep(1, K), list(), list())
  shapes2_t <- list(rep(1, K), list(), list())
  shapes1_a <- rep(1, N)
  shapes2_a <- rep(1, N)
  
  for (i in 1:K){
    idx <- which(cluster_level1 == i)
    cluster_level2[[i]] <- kmeans(xs[idx], K)$cluster
    cluster_level3[[i]] <- list()
    tree[[i]] <- list()
    shapes1_t[[2]][[i]] <- rep(1, K)
    shapes2_t[[2]][[i]] <- rep(1, K)
    shapes1_t[[3]][[i]] <- list()
    shapes2_t[[3]][[i]] <- list()
    
    for (j in 1:K){
      idx2 <- which(cluster_level2[[i]] == j)
      tree[[i]][[j]] <- list()
      
      if (length(idx2) > K){
        cluster_level3[[i]][[j]] <- kmeans(xs[idx][idx2], K)$cluster
        shapes1_t[[3]][[i]][[j]] <- rep(1, K)
        shapes2_t[[3]][[i]][[j]] <- rep(1, K)
        for (k in 1:K){
          tree[[i]][[j]][[k]] <- idx[idx2][which(cluster_level3[[i]][[j]] == k)]
        }
      } else {
        cluster_level3[[i]][[j]] <- 1:length(idx2)
        shapes1_t[[3]][[i]][[j]] <- rep(1, length(idx2))
        shapes2_t[[3]][[i]][[j]] <- rep(1, length(idx2))
        for (k in 1:length(idx2)){
          tree[[i]][[j]][[k]] <- idx[idx2][which(cluster_level3[[i]][[j]] == k)]
        }
      }
    }
  }
  
  rm(cluster_level1, cluster_level2, cluster_level3)
  
  shapes1_t <- shapes1_t[1:L]
  shapes2_t <- shapes2_t[1:L]
  
  if (L <= 2){
    for (i in 1:length(tree)){
      for (j in 1:length(tree[[i]])){
        tree[[i]][[j]] <- unlist(tree[[i]][[j]])
      }
    }
    if (L <= 1){
      for (i in 1:length(tree)){
        tree[[i]] <- unlist(tree[[i]])
      }
    }
  }
  
  return(list(shapes1_a=shapes1_a, shapes2_a=shapes2_a,
              shapes1_t=shapes1_t, shapes2_t=shapes2_t,
              tree=tree, theta=theta))
}

hts_non_binary <- function(shapes1_t, shapes2_t, shapes1_a, shapes2_a, theta, tree, N_max=3000){
  # Hierarchical Thompson sampling algorithm for non-binary trees
  
  rewards <- c()
  L <- length(shapes1_t)
  
  for (i in 1:N_max){
    samp <- rbeta(length(tree), shapes1_t[[1]], shapes2_t[[1]])
    sel_node <- order(samp, decreasing=TRUE)[1]
    path <- c(sel_node)
    node <- tree[[sel_node]]
    
    if (L >= 2){
      samp <- rbeta(length(node), shapes1_t[[2]][[sel_node]], shapes2_t[[2]][[sel_node]])
      sel_node <- order(samp, decreasing=TRUE)[1]
      path <- append(path, sel_node)
      node <- node[[sel_node]]
      
      if (L >= 3){
        samp <- rbeta(length(node), shapes1_t[[3]][[path[1]]][[sel_node]], shapes2_t[[3]][[path[1]]][[sel_node]])
        sel_node <- order(samp, decreasing=TRUE)[1]
        path <- append(path, sel_node)
        node <- node[[sel_node]]
      }
    }
    
    samp <- rbeta(length(node), shapes1_a[node], shapes2_a[node])
    a <- node[order(samp, decreasing=TRUE)[1]]
    r <- rbinom(1, 1, theta[a])
    rewards[i] <- r
    
    shapes1_a[a] <- shapes1_a[a] + r
    shapes2_a[a] <- shapes2_a[a] + (1 - r)
    shapes1_t[[1]][path[1]] <- shapes1_t[[1]][path[1]] + r
    shapes2_t[[1]][path[1]] <- shapes2_t[[1]][path[1]] + (1 - r)
    if (L >= 2){
      shapes1_t[[2]][[path[1]]][path[2]] <- shapes1_t[[2]][[path[1]]][path[2]] + r
      shapes2_t[[2]][[path[1]]][path[2]] <- shapes2_t[[2]][[path[1]]][path[2]] + (1 - r)
      
      if (L >= 3){
        shapes1_t[[3]][[path[1]]][[path[2]]][path[3]] <- shapes1_t[[3]][[path[1]]][[path[2]]][path[3]] + r
        shapes2_t[[3]][[path[1]]][[path[2]]][path[3]] <- shapes2_t[[3]][[path[1]]][[path[2]]][path[3]] + (1 - r)
      }
    }
  }
  
  regret <- max(theta) - rewards
  return(cumsum(regret))
}

uct <- function(played_node, played_arm, theta, tree, N_max=3000){
  # Upper cofidence bound algorithm for trees (Kocsis et al. 2006)
  # (UCB1 algorithm applied at each node)
  # Currently implemented for max 3 levels L in the hierarchical structure
  
  rewards <- c()
  node_r <- played_node
  arm_r <- played_arm
  L <- length(played_node)
  
  for (i in 1:N_max){
    samp <- (node_r[[1]] - 1) / played_node[[1]] + sqrt(2 * log(i) / played_node[[1]]) # step into node that maximizes this expression
    sel_node <- order(samp, decreasing=TRUE)[1]
    path <- c(sel_node)
    node <- tree[[sel_node]]
    
    if (L >= 2){
      samp <- (node_r[[2]][[sel_node]] - 1) / played_node[[2]][[sel_node]] + sqrt(2 * log(played_node[[1]][sel_node]) / played_node[[2]][[sel_node]])
      sel_node <- order(samp, decreasing=TRUE)[1]
      path <- append(path, sel_node)
      node <- node[[sel_node]]
      
      if (L >= 3){
        samp <- (node_r[[3]][[path[1]]][[sel_node]] - 1) / played_node[[3]][[path[1]]][[sel_node]] + sqrt(2 * log(played_node[[2]][[path[1]]][sel_node]) / played_node[[3]][[path[1]]][[sel_node]])
        sel_node <- order(samp, decreasing=TRUE)[1]
        path <- append(path, sel_node)
        node <- node[[sel_node]]
      }
    }
    
    samp <- (arm_r[node] - 1) / played_arm[node] + sqrt(2 * log(sum(played_arm[node])) / played_arm[node])
    a <- node[order(samp, decreasing=TRUE)[1]]
    r <- rbinom(1, 1, theta[a])
    rewards[i] <- r
    
    played_arm[a] <- played_arm[a] + 1
    arm_r[a] <- arm_r[a] + r
    played_node[[1]][path[1]] <- played_node[[1]][path[1]] + 1
    node_r[[1]][path[1]] <- node_r[[1]][path[1]] + r
    if (L >= 2){
      played_node[[2]][[path[1]]][path[2]] <- played_node[[2]][[path[1]]][path[2]] + 1
      node_r[[2]][[path[1]]][path[2]] <- node_r[[2]][[path[1]]][path[2]] + r
      
      if (L >= 3){
        played_node[[3]][[path[1]]][[path[2]]][path[3]] <- played_node[[3]][[path[1]]][[path[2]]][path[3]] + 1
        node_r[[3]][[path[1]]][[path[2]]][path[3]] <- node_r[[3]][[path[1]]][[path[2]]][path[3]] + r
      }
    }
  }
  
  regret <- max(theta) - rewards
  return(cumsum(regret))
}

### Clustered contextual bandits

prepare_contextual_cluster <- function(N, K, eps, d=5, seed=1, random=FALSE){
  # Assign clusters and probabilities to each of the arms
  # Prepare data to be processed with Thompson Sampling algorithms
  
  set.seed(seed)
  
  clusters <- sample(1:K, N, replace=TRUE)
  
  theta_c <- rmvnorm(K, mean=rep(0, d), sigma=diag(1, d)) # cluster centroids
  v_a <- rmvnorm(N, mean=rep(0, d), sigma=diag(1, d)) # arm difference to its centroid
  
  theta_a <- theta_c[clusters,] + eps * v_a # arm coefficients
  
  if (random) {
    # Randomly resample assigned clusters - used in results beyond original paper section
    clusters <- sample(1:K, N, replace=TRUE)
  }
  
  return(list(theta_a=theta_a, clusters=clusters))
}

lin_ts <- function(theta_a, N_max=1000){
  # LinTS algorithm for Thompson sampling in linear contextual bandits (Agrawal 2012)
  
  regret <- c()
  N <- dim(theta_a)[1]
  d <- dim(theta_a)[2]
  B_a <- list()
  f_a <- matrix(0, nrow=N, ncol=d)
  mu_a <- matrix(0, nrow=N, ncol=d)
  
  for (i in 1:N){
    B_a[[i]] <- diag(1, d)
  }
  
  for (i in 1:N_max){
    x <- rmvnorm(1, mean=rep(0, d), sigma=diag(1, d))# context
    opt_rew <- max(theta_a %*% t(x))
    mu <- mu_a %*% t(x) # means
    sigma <- c()
    for (j in 1:N){
      sigma[j] <- x %*% solve(B_a[[j]]) %*% t(x) #standard deviation
    }
    samp_a <- rnorm(N, mu, sigma)
    a <- order(samp_a, decreasing=TRUE)[1]
    r <- sign(theta_a[a,] %*% t(x)) * runif(1, 0, abs(2 * theta_a[a,] %*% t(x)))
    regret[i] <- opt_rew - r
    B_a[[a]] <- B_a[[a]] + t(x) %*% x
    f_a[a,] <- f_a[a,] + r %*% x
    mu_a[a,] <- solve(B_a[[a]]) %*% f_a[a,]
  }
  return(cumsum(regret))
}

lin_tsc <- function(clusters, theta_a, N_max=1000){
  # LinTSC algorithm for Thompson sampling in linear contextual bandits with clustered arms
  
  regret <- c()
  N <- dim(theta_a)[1]
  d <- dim(theta_a)[2]
  K <- length(unique(clusters))
  B_a <- list()
  B_c <- list()
  f_a <- matrix(0, nrow=N, ncol=d)
  f_c <- matrix(0, nrow=K, ncol=d)
  mu_a <- matrix(0, nrow=N, ncol=d)
  mu_c <- matrix(0, nrow=K, ncol=d)
  
  for (i in 1:K){
    B_c[[i]] <- diag(1, d)
  }
  
  for (i in 1:N){
    B_a[[i]] <- diag(1, d)
  }
  
  for (i in 1:N_max){
    x <- rmvnorm(1, mean=rep(0, d), sigma=diag(1, d)) # context
    opt_rew <- max(theta_a %*% t(x))
    mu_cluster <- mu_c %*% t(x) # cluster means
    sigma_cluster <- c()
    for (j in 1:K){
      sigma_cluster[j] <- x %*% solve(B_c[[j]]) %*% t(x) # cluster standard deviation
    }
    samp_c <- rnorm(K, mu_cluster, sigma_cluster)
    k <- order(samp_c, decreasing=TRUE)[1]
    idx <- which(clusters == k)
    mu <- mu_a[idx,] %*% t(x) # arm means
    sigma <- c()
    for (j in 1:length(idx)){
      sigma[j] <- x %*% solve(B_a[[idx[j]]]) %*% t(x) # arm standard deviation
    }
    samp_a <- rnorm(length(idx), mu, sigma)
    a <- idx[order(samp_a, decreasing=TRUE)[1]]
    r <- sign(theta_a[a,] %*% t(x)) * runif(1, 0, abs(2 * theta_a[a,] %*% t(x)))
    regret[i] <- opt_rew - r
    B_a[[a]] <- B_a[[a]] + t(x) %*% x
    B_c[[k]] <- B_c[[k]] + t(x) %*% x
    f_a[a,] <- f_a[a,] + r %*% x
    f_c[k,] <- f_c[k,] + r %*% x
    mu_a[a,] <- solve(B_a[[a]]) %*% f_a[a,]
    mu_c[k,] <- solve(B_c[[k]]) %*% f_c[k,]
  }
  return(cumsum(regret))
}

lin_ucb <- function(theta_a, N_max=3000){
  # LinUCB algorithm for contextual bandits (Li et al. 2010)
  
  regret <- c()
  N <- dim(theta_a)[1]
  d <- dim(theta_a)[2]
  B_a <- list()
  f_a <- matrix(0, nrow=N, ncol=d)
  mu_a <- matrix(0, nrow=N, ncol=d)
  
  for (i in 1:N){
    B_a[[i]] <- diag(1, d)
  }
  
  for (i in 1:N_max){
    x <- rmvnorm(1, mean=rep(0, d), sigma=diag(1, d)) # context
    opt_rew <- max(theta_a %*% t(x))
    p <- mu_a %*% t(x)
    for (j in 1:N){
      p[j] <- p[j] + 2 * x %*% solve(B_a[[j]]) %*% t(x) # upper confidence bound vector
    }
    a <- order(p, decreasing=TRUE)[1]
    r <- sign(theta_a[a,] %*% t(x)) * runif(1, 0, abs(2 * theta_a[a,] %*% t(x)))
    regret[i] <- opt_rew - r
    B_a[[a]] <- B_a[[a]] + t(x) %*% x
    f_a[a,] <- f_a[a,] + r %*% x
    mu_a[a,] <- solve(B_a[[a]]) %*% f_a[a,]
  }
  return(cumsum(regret))
}

lin_ucbc <- function(clusters, theta_a, N_max=1000){
  # LinUCBC algorithm for clustered contextual bandits (Bouneffouf et al. 2019)
  
  regret <- c()
  N <- dim(theta_a)[1]
  d <- dim(theta_a)[2]
  K <- length(unique(clusters))
  B_a <- list()
  B_c <- list()
  f_a <- matrix(0, nrow=N, ncol=d)
  f_c <- matrix(0, nrow=K, ncol=d)
  mu_a <- matrix(0, nrow=N, ncol=d)
  mu_c <- matrix(0, nrow=K, ncol=d)
  
  for (i in 1:K){
    B_c[[i]] <- diag(1, d)
  }
  
  for (i in 1:N){
    B_a[[i]] <- diag(1, d)
  }
  
  for (i in 1:N_max){
    x <- rmvnorm(1, mean=rep(0, d), sigma=diag(1, d)) # context
    opt_rew <- max(theta_a %*% t(x))
    p_c <- mu_c %*% t(x)
    for (j in 1:K){
      p_c[j] <- p_c[j] + 2 * x %*% solve(B_c[[j]]) %*% t(x) # upper confidence bound vector - clusters
    }
    k <- order(p_c, decreasing=TRUE)[1] # select most promising cluster
    idx <- which(clusters == k)
    p_a <- mu_a[idx,] %*% t(x)
    for (j in 1:length(idx)){
      p_a[j] <- p_a[j] + 2 * x %*% solve(B_a[[idx[j]]]) %*% t(x) # upper confidence bound vector - arms
    }
    a <- idx[order(p_a, decreasing=TRUE)[1]] # select most promising arm
    r <- sign(theta_a[a,] %*% t(x)) * runif(1, 0, abs(2 * theta_a[a,] %*% t(x)))
    regret[i] <- opt_rew - r
    B_a[[a]] <- B_a[[a]] + t(x) %*% x
    B_c[[k]] <- B_c[[k]] + t(x) %*% x
    f_a[a,] <- f_a[a,] + r %*% x
    f_c[k,] <- f_c[k,] + r %*% x
    mu_a[a,] <- solve(B_a[[a]]) %*% f_a[a,]
    mu_c[k,] <- solve(B_c[[k]]) %*% f_c[k,]
  }
  return(cumsum(regret))
}