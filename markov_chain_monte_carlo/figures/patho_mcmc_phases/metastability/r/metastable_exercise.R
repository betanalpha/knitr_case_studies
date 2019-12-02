############################################################
#
# Initial setup
#
############################################################

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_trans <- c("#DCBCBC30")
c_dark_trans <- c("#8F272730")

library(rstan)

fast_monitor <- function(sims, warmup, names=NULL) {
  dim_sims <- dim(sims)
  if (is.null(dim_sims)) {
    dim(sims) <- c(length(sims), 1, 1)
  } else if (length(dim_sims) == 2) {
    dim(sims) <- c(dim_sims, 1)
  } else if (length(dim_sims) > 3) {
    stop("'sims' has more than 3 dimensions")
  }
  
  parnames <- names
  if (is.null(parnames)) {
    parnames <- paste0("V", seq_len(dim(sims)[3]))
  }
  iter <- dim(sims)[1]
  chains <- dim(sims)[2]
  if (warmup > dim(sims)[1]) {
    stop("warmup is larger than the total number of iterations")
  }
  if (warmup >= 1) {
    sims <- sims[-seq_len(warmup), , , drop = FALSE]
  }
  
  out <- vector("list", length(parnames))
  out <- setNames(out, parnames)
  
  for (i in seq_along(out)) {
    sims_i <- sims[, , i]
    mean <- mean(sims_i)
    var <- var(as.vector(sims_i))
    mcse_mean <- rstan:::mcse_mean(sims_i)
    ess <- round(rstan:::ess_rfun(sims_i))
    rhat <- rstan:::rhat_rfun(rstan:::split_chains(sims_i))
    out[[i]] <- c(mean, var, mcse_mean, ess, rhat)
  }
  
  out <- as.data.frame(do.call(rbind, out))
  colnames(out) <- c("mean", "var", "se_mean", "n_eff", "Rhat")
  rownames(out) <- parnames
  
  out <- structure(out, chains = chains, iter = iter, warmup = warmup,
                   class = c("simsummary", "data.frame"))
  return(invisible(out))
}

############################################################
#
# How does the Random Walk Metropolis algorithm perform
# on a target distribution with a two-dimensional Gaussian 
# density function?
#
############################################################

D <- 2

target_lpdf1 <- function(q) {
  mu.1 <- 4
  sigma.1 <- 1
  
  mu.2 <- 8
  sigma.2 <- 2
  
  lpdf <- -0.5 * ( ((q[1] - mu.1) / sigma.1)**2 + ((q[2] - mu.2) / sigma.2)**2 ) 
  lpdf <- lpdf - 0.5 * log(6.283185307179586) - log(sigma.1) - log(sigma.2)
  lpdf
}

target_lpdf2 <- function(q) {
  mu.1 <- -8
  sigma.1 <- 2
  
  mu.2 <- -4
  sigma.2 <- 1
  
  lpdf <- -0.5 * ( ((q[1] - mu.1) / sigma.1)**2 + ((q[2] - mu.2) / sigma.2)**2 ) 
  lpdf <- lpdf - 0.5 * log(6.283185307179586) - log(sigma.1) - log(sigma.2)
  lpdf
}

target_lpdf <- function(x) {
  lpdf1 <- log(0.5) + target_lpdf1(x)
  lpdf2 <- log(0.5) + target_lpdf2(x)
  if (lpdf1 > lpdf2) {
    lpdf <- lpdf1 + log(1 + exp(lpdf2 - lpdf1))
  } else {
    lpdf <- lpdf2 + log(1 + exp(lpdf1 - lpdf2))
  }
  (lpdf)
}

N <- 200
q.1 <- seq(-13, 7, 20 / N)
q.2 <- seq(-7, 13, 20 / N)
densities <- matrix(data = 0, nrow = N + 1, ncol = N + 1)

for (n in 1:N) {
  for (m in 1:N) {
    q <- c(q.1[n], q.2[m])
    densities[n, m] <- exp(target_lpdf(q)) 
  }
}

contour(q.1, q.2, densities, nlevels=30, main="Multimodal Target Density",
        xlab="q.1", xlim=c(-13, 7), 
        ylab="q.2", ylim=c(-7, 13), drawlabels = FALSE, col = c_dark, lwd = 2)

# Fit
n_transitions <- 5000
set.seed(5849586)

sigma <- 2

mcmc_samples <- matrix(data = 0, nrow = n_transitions + 1, ncol = D + 1)

# Seed the Markov chain from a diffuse sample
mcmc_samples[1, 1:D] <- rnorm(D, 0, 5) 
mcmc_samples[1, D + 1] <- 1

for (n in 1:n_transitions) {
  q0 <- mcmc_samples[n, 1:D] # Initial point
  qp <- rnorm(D, q0, sigma)  # Proposal
  
  # Compute acceptance probability
  accept_prob <- min(1, exp(target_lpdf(qp) - target_lpdf(q0)))
  mcmc_samples[n, D + 1] <- accept_prob
  
  # Apply Metropolis correction
  u = runif(1, 0, 1)
  if (accept_prob > u)
    mcmc_samples[n + 1, 1:D] <- qp
  else
    mcmc_samples[n + 1, 1:D] <- q0
}

parameter_names <- c("varpi_1", "varpi_2", "accept_prob")
mcmc_stats <- fast_monitor(array(mcmc_samples, dim = c(n_transitions + 1, 1, 3)), 
                           warmup=100, parameter_names)
mcmc_stats

par(mfrow=c(1, 2))

plot(0:n_transitions, mcmc_samples[,1], main="",
     type="l", lwd=0.5, col=c_dark,
     xlab="Iteration", xlim=c(0, n_transitions),
     ylab="varpi_1")

plot(0:n_transitions, mcmc_samples[,2], main="",
     type="l", lwd=0.5, col=c_light,
     xlab="Iteration", xlim=c(0, n_transitions),
     ylab="varpi_2")
  
# One long chain
n_transitions <- 200000
set.seed(5849586)

mcmc_samples <- matrix(data = 0, nrow = n_transitions + 1, ncol = D + 1)

mcmc_samples[1, 1:D] <- rnorm(D, 0, 5)
mcmc_samples[1, D + 1] <- 1

for (n in 1:n_transitions) {
  q0 <- mcmc_samples[n, 1:D]
  qp <- rnorm(D, q0, sigma)
  
  accept_prob <- min(1, exp(target_lpdf(qp) - target_lpdf(q0)))
  mcmc_samples[n + 1, D + 1] <- accept_prob
  
  u = runif(1, 0, 1)
  if (accept_prob > u)
    mcmc_samples[n + 1, 1:D] <- qp
  else
    mcmc_samples[n + 1, 1:D] <- q0
}

mcmc_stats <- fast_monitor(array(mcmc_samples, dim = c(n_transitions + 1, 1, D + 1)), 
                           warmup=100, parameter_names)
mcmc_stats

par(mfrow=c(1, 2))

plot(0:n_transitions, mcmc_samples[,1], main="",
     type="l", lwd=0.5, col=c_dark,
     xlab="Iteration", xlim=c(0, n_transitions),
     ylab="varpi_1")

plot(0:n_transitions, mcmc_samples[,2], main="",
     type="l", lwd=0.5, col=c_light,
     xlab="Iteration", xlim=c(0, n_transitions),
     ylab="varpi_2")

# Many chains
n_chains <- 4
n_transitions <- 5000
set.seed(5849586)

mcmc_samples = array(0, dim=c(n_transitions + 1, 4, D + 1))

for (c in 1:n_chains) {
  mcmc_samples[1, c, 1:D] <- rnorm(D, 0, 5)
  mcmc_samples[1, c, D + 1] <- 1
  
  for (n in 1:n_transitions) {
    q0 <- mcmc_samples[n, c, 1:D]
    qp <- rnorm(D, q0, sigma)
    
    accept_prob <- min(1, exp(target_lpdf(qp) - target_lpdf(q0)))
    mcmc_samples[n + 1, c, D + 1] <- accept_prob
    
    u = runif(1, 0, 1)
    if (accept_prob > u)
      mcmc_samples[n + 1, c, 1:D] <- qp
    else
      mcmc_samples[n + 1, c, 1:D] <- q0
  }
}

mcmc_stats <- fast_monitor(mcmc_samples, warmup=100, parameter_names)
mcmc_stats

par(mfrow=c(2, 4))
for (c in 1:4) {
  plot(0:n_transitions, mcmc_samples[,c,1], main=paste("Chain", c, sep=" "),
       type="l", lwd=0.5, col=c_dark,
       xlab="Iteration", xlim=c(0, n_transitions),
       ylab="varpi_1", ylim=c(-13, 7))
}

for (c in 1:4) {  
  plot(0:n_transitions, mcmc_samples[,c,2], main=paste("Chain", c, sep=" "),
       type="l", lwd=0.5, col=c_light,
       xlab="Iteration", xlim=c(0, n_transitions),
       ylab="varpi_2", c(-7, 13))
}


par(mfrow=c(1, 1))
contour(q.1, q.2, densities, nlevels=30, main="Metastable Target Density",
        xlab="q.1", xlim=c(-13, 7), 
        ylab="q.2", ylim=c(-7, 13), drawlabels = FALSE, col = c_dark, lwd = 2)
lines(mcmc_samples[,1,1], mcmc_samples[,1,2], lwd=0.5, col=c_light)