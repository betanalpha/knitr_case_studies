############################################################
# Initial setup
############################################################

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

util <- new.env()
source('stan_utility.R', local=util)
source('plot_utility.R', local=util)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_trans <- c("#DCBCBC80")
c_dark_trans <- c("#8F272780")
c_green_trans <- c("#00FF0080")

par(family="CMU Serif", las=1, bty="l", cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 1))

set.seed(58533858)
#set.seed(78533858)

############################################################
#
# Simulate Data
#
############################################################

# Let's simulate some data with a varying distribution of parameter sizes.
# Three are around the nominal scale ot 1, three are smaller, and three 
# are larger.

K <- 12

theta_true <- c(-0.09604655,  0.02063505,  0.01246716, 
                -0.04487767, -0.13519031,  0.09421304,
                0.061588542,  0.004749712, -0.177075605,
                -29.449233,  32.997172,  18.517443)
N <- K
sigma <- 1
context_idx <- rep(1:K, 1)

y <- rnorm(K, theta_true[context_idx], sigma)

data <- list("N" = N, "K" = K, "context_idx" = context_idx,
             "y" = y, "sigma" = sigma)

############################################################
#
# Normal
#
############################################################

# Narrow
fit <- stan(file='stan_programs/normal_narrow.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

samples = extract(fit)

par(mfrow=c(3, 3))

for (k in 1:9) {
  hist(samples$theta[, k], breaks=seq(-40, 40, 0.25),
       main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
       xlab="theta", yaxt='n', ylab="")
  abline(v=theta_true[k], col="white", lwd=2)
  abline(v=theta_true[k], col="black", lwd=1)
}

# Wide
fit <- stan(file='stan_programs/normal_wide.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

samples = extract(fit)

par(mfrow=c(3, 3))

for (k in 1:9) {
  hist(samples$theta[, k], breaks=seq(-40, 40, 0.25),
       main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
       xlab="theta", yaxt='n', ylab="")
  abline(v=theta_true[k], col="white", lwd=2)
  abline(v=theta_true[k], col="black", lwd=1)
}

par(mfrow=c(1, 1))

k <- 1
hist(samples$theta[, k], breaks=seq(-5, 5, 0.15),
     main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
     xlab="theta", yaxt='n', ylab="")
abline(v=theta_true[k], col="white", lwd=2)
abline(v=theta_true[k], col="black", lwd=1)

# Hierarchical
fit <- stan(file='stan_programs/hier_normal_cp.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

samples = extract(fit)

par(mfrow=c(1, 1))

hist(abs(rnorm(4000, 0, 10)), breaks=seq(0, 50, 0.5),
     main="", col=c_light, border=c_light_highlight,
     xlab="tau", yaxt='n', ylab="", ylim=c(0, 300))
hist(samples$tau, breaks=seq(0, 50, 0.5),
     main="", col=c_dark, border=c_dark_highlight, add=T)

abline(v=0.1, col="white", lwd=2)
abline(v=0.1, col="black", lwd=1)

abline(v=20, col="white", lwd=2)
abline(v=20, col="black", lwd=1)

par(mfrow=c(3, 3))

for (k in 1:9) {
  hist(samples$theta[, k], breaks=seq(-40, 40, 0.25),
       main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
       xlab="theta", yaxt='n', ylab="")
  abline(v=theta_true[k], col="white", lwd=2)
  abline(v=theta_true[k], col="black", lwd=1)
}

par(mfrow=c(1, 1))

k <- 1
hist(samples$theta[, k], breaks=seq(-5, 5, 0.15),
     main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
     xlab="theta", yaxt='n', ylab="")
abline(v=theta_true[k], col="white", lwd=2)
abline(v=theta_true[k], col="black", lwd=1)

############################################################
#
# Normal Mixture
#
############################################################

seed1 <- .Random.seed
.Random.seed <- seed1

# Full centered parameterization
fit <- stan(file='stan_programs/hier_normal_mixture_cp.stan', data=data, 
            seed=4938483, refresh=1000)

simu_inits <- function(chain_id) {
  return(list("tau1" = abs(rnorm(1, 0, 0.1))))
}

fit <- stan(file='stan_programs/hier_normal_mixture_cp.stan', data=data, 
            seed=4938483, refresh=1000, init=simu_inits)

util$check_all_diagnostics(fit)

sapply(1:4, function(c) get_sampler_params(fit, inc_warmup=FALSE)[[c]][,'stepsize__'][1])

sapply(1:4, function(c) range(get_sampler_params(fit, inc_warmup=FALSE)[[c]][,'treedepth__']))


partition <- util$partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mfrow=c(3, 3))

for (k in 1:9) {
  name_x <- paste("theta[", k, "]", sep='')
  name_y <- "tau1"
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("k = ", k),
       xlab=name_x, xlim=c(theta_true[k] - 5, theta_true[k] + 5), 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-7, 0))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
}

par(mfrow=c(3, 3))

for (k in 1:9) {
  name_x <- paste("theta[", k, "]", sep='')
  name_y <- "tau2"
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("k = ", k),
       xlab=name_x, xlim=c(0 - 5, 0 + 5), 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-7, 4.5))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
}

unpermuted_samples <- extract(fit, permute=FALSE)

par(mfrow=c(2, 2))

plot_idx <- 7

for (c in 1:4) {
  plot(1:1000, unpermuted_samples[,c,plot_idx], type="l", lwd=1, col=c_dark,
       main=paste("Chain", c),
       xlab="Iteration",  xlim=c(1, 1000),
       ylab="thetea_ncp[4]", ylim=c(-500, 500))
}

par(mfrow=c(2, 2))

plot_idx <- data$K + 2

for (c in 1:4) {
  plot(1:1000, log(unpermuted_samples[,c,plot_idx]), type="l", lwd=1, col=c_dark,
       main=paste("Chain", c),
       xlab="Iteration",  xlim=c(1, 1000),
       ylab="tau1", ylim=c(-10, 0))
}

par(mfrow=c(2, 2))

plot_idx <- 12

for (c in 1:4) {
  plot(1:1000, log(unpermuted_samples[,c,plot_idx]), type="l", lwd=1, col=c_dark,
       main=paste("Chain", c),
       xlab="Iteration",  xlim=c(1, 1000),
       ylab="tau2", ylim=c(-7, 4))
}

# Partially center stuff?
data$w <- c(0, 0, 0,
            0, 0, 0,
            1, 1, 1)

fit <- stan(file='stan_programs/hier_normal_mixture_partial.stan', 
            data=data, seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

partition <- util$partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

# Look at the overlay of the funnels!

par(mfrow=c(3, 3))

for (k in 1:9) {
  name_x <- paste("theta_tilde[", k, "]", sep='')
  name_y <- "tau1"
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("k = ", k),
       xlab=name_x, xlim=c(theta_true[k] - 50, theta_true[k] + 50), 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-9, 0))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
}

sapply(1:4, function(c) get_sampler_params(fit, inc_warmup=FALSE)[[c]][,'stepsize__'][1])

sapply(1:4, function(c) range(get_sampler_params(fit, inc_warmup=FALSE)[[c]][,'treedepth__']))

data$w <- c(0.2, 0.2, 0.2,
            0.2, 0.2, 0.2,
            1, 1, 1)

fit <- stan(file='stan_programs/hier_normal_mixture_partial.stan', 
            data=data, seed=4938483, refresh=1000, init_r=1)

util$check_n_eff(fit)
util$check_rhat(fit)
util$check_div(fit)
util$check_treedepth(fit, 15)
util$check_energy(fit)

sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
range(treedepths)
which(treedepths==15)

# Non-center the smaller parameters
seed2 <- .Random.seed

data$cp_idx <- c(7, 8, 9)
data$K_cp <- length(data$cp_idx)

data$ncp_idx <- setdiff(1:data$K, data$cp_idx)
data$K_ncp <- length(data$ncp_idx)

fit <- stan(file='stan_programs/hier_normal_mixture_mixed.stan', data=data, 
            seed=4938483, refresh=1000, init=simu_inits)

util$check_all_diagnostics(fit)

partition <- util$partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mfrow=c(3, 3))

for (k in 1:9) {
  if (k %in% data$cp_idx)
    name_x <- paste("theta_cp[", which(data$cp_idx == k), "]", sep='')
  else
    name_x <- paste("theta_ncp[", which(data$ncp_idx == k), "]", sep='')
  
  name_y <- "tau1"
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("k = ", k),
       xlab=name_x, xlim=c(theta_true[k] - 50, theta_true[k] + 50), 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-9, 0))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
}


par(mfrow=c(3, 3))

for (k in 1:9) {
  if (k %in% data$cp_idx)
    name_x <- paste("theta_cp[", which(data$cp_idx == k), "]", sep='')
  else
    name_x <- paste("theta_ncp[", which(data$ncp_idx == k), "]", sep='')
  
  name_y <- "gamma"
  
  plot(nondiv_params[name_x][,1], nondiv_params[name_y][,1],
       col=c_dark_trans, pch=16, main=paste("k = ", k),
       xlab=name_x, xlim=c(theta_true[k] - 50, theta_true[k] + 50), 
       ylab=paste("", name_y, "", sep=""), ylim=c(0, 1))
  points(div_params[name_x][,1], div_params[name_y][,1],
         col=c_green_trans, pch=16)
}


par(mfrow=c(3, 3))

for (k in 1:9) {
  if (k %in% data$cp_idx)
    name_x <- paste("theta_cp[", which(data$cp_idx == k), "]", sep='')
  else
    name_x <- paste("theta_ncp[", which(data$ncp_idx == k), "]", sep='')
  
  name_y <- "tau2"
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("k = ", k),
       xlab=name_x, xlim=c(0 - 100, 0 + 100), 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-4, 4))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
}

name_x <- "gamma"
name_y <- "tau1"

plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
     col=c_dark_trans, pch=16, main=paste("k = ", k),
     xlab=name_x, xlim=c(0, 1), 
     ylab=paste("log(", name_y, ")", sep=""), ylim=c(-8, 4))
points(div_params[name_x][,1], log(div_params[name_y][,1]),
       col=c_green_trans, pch=16)

name_x <- "gamma"
name_y <- "tau2"

plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
     col=c_dark_trans, pch=16, main=paste("k = ", k),
     xlab=name_x, xlim=c(0, 1), 
     ylab=paste("log(", name_y, ")", sep=""), ylim=c(-4, 4))
points(div_params[name_x][,1], log(div_params[name_y][,1]),
       col=c_green_trans, pch=16)

samples = extract(fit)

par(mfrow=c(1, 2))

hist(abs(rnorm(4000, 0, 0.1)), breaks=seq(0, 0.5, 0.01),
     main="", col=c_light, border=c_light_highlight,
     xlab="tau1", yaxt='n', ylab="", ylim=c(0, 400))
hist(samples$tau1, breaks=seq(0, 0.5, 0.01),
     main="", col=c_dark, border=c_dark_highlight, add=T)

abline(v=0.1, col="white", lwd=2)
abline(v=0.1, col="black", lwd=1)

hist(abs(rnorm(4000, 0, 10)), breaks=seq(0, 50, 1),
     main="", col=c_light, border=c_light_highlight,
     xlab="tau2", yaxt='n', ylab="", ylim=c(0, 400))
hist(samples$tau2, breaks=seq(0, 50, 1),
     main="", col=c_dark, border=c_dark_highlight, add=T)

abline(v=20, col="white", lwd=2)
abline(v=20, col="black", lwd=1)

par(mfrow=c(9, 2))

for (k in 1:9) {
  hist(samples$theta[, k], breaks=seq(theta_true[k] - 5, theta_true[k] + 5, 0.1),
       main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
       xlab="theta", yaxt='n', ylab="")
  abline(v=theta_true[k], col="white", lwd=2)
  abline(v=theta_true[k], col="black", lwd=1)
  
  hist(samples$lambda[, k], breaks=seq(0, 1, 0.01),
       main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
       xlab="lambda", yaxt='n', ylab="")
}

############################################################
#
# Laplace
#
############################################################

# Narrow
fit <- stan(file='laplace_narrow.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

samples = extract(fit)

par(mfrow=c(3, 3))

for (k in 1:9) {
  hist(samples$theta[, k], breaks=seq(-35, 35, 0.25),
       main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
       xlab="theta", yaxt='n', ylab="")
  abline(v=theta_true[k], col="white", lwd=2)
  abline(v=theta_true[k], col="black", lwd=1)
}

# Wide
fit <- stan(file='laplace_wide.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

samples = extract(fit)

par(mfrow=c(3, 3))

for (k in 1:9) {
  hist(samples$theta[, k], breaks=seq(-40, 40, 0.25),
       main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
       xlab="theta", yaxt='n', ylab="")
  abline(v=theta_true[k], col="white", lwd=2)
  abline(v=theta_true[k], col="black", lwd=1)
}

par(mfrow=c(1, 1))

k <- 1
hist(samples$theta[, k], breaks=seq(-5, 5, 0.15),
     main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
     xlab="theta", yaxt='n', ylab="")
abline(v=theta_true[k], col="white", lwd=2)
abline(v=theta_true[k], col="black", lwd=1)

# Hierarchical
fit <- stan(file='hier_laplace_cp.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

samples = extract(fit)

par(mfrow=c(1, 1))

hist(abs(rnorm(4000, 0, 10)), breaks=seq(0, 50, 0.5),
     main="", col=c_light, border=c_light_highlight,
     xlab="tau", yaxt='n', ylab="", ylim=c(0, 300))
hist(samples$tau, breaks=seq(0, 50, 0.5),
     main="", col=c_dark, border=c_dark_highlight, add=T)

abline(v=0.1, col="white", lwd=2)
abline(v=0.1, col="black", lwd=1)

abline(v=20, col="white", lwd=2)
abline(v=20, col="black", lwd=1)

par(mfrow=c(3, 3))

for (k in 1:9) {
  hist(samples$theta[, k], breaks=seq(-40, 40, 0.25),
       main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
       xlab="theta", yaxt='n', ylab="")
  abline(v=theta_true[k], col="white", lwd=2)
  abline(v=theta_true[k], col="black", lwd=1)
}

par(mfrow=c(1, 1))

k <- 1
hist(samples$theta[, k], breaks=seq(-5, 5, 0.15),
     main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
     xlab="theta", yaxt='n', ylab="")
abline(v=theta_true[k], col="white", lwd=2)
abline(v=theta_true[k], col="black", lwd=1)

############################################################
#
# Cauchy
#
############################################################

# Narrow
fit <- stan(file='cauchy_narrow.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

samples = extract(fit)

par(mfrow=c(3, 3))

for (k in 1:9) {
  hist(samples$theta[, k], breaks=seq(-40, 40, 0.25),
       main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
       xlab="theta", yaxt='n', ylab="")
  abline(v=theta_true[k], col="white", lwd=2)
  abline(v=theta_true[k], col="black", lwd=1)
}

# Wide
fit <- stan(file='cauchy_wide.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

samples = extract(fit)

par(mfrow=c(3, 3))

for (k in 1:9) {
  hist(samples$theta[, k], breaks=seq(-40, 40, 0.25),
       main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
       xlab="theta", yaxt='n', ylab="")
  abline(v=theta_true[k], col="white", lwd=2)
  abline(v=theta_true[k], col="black", lwd=1)
}

par(mfrow=c(1, 1))

k <- 1
hist(samples$theta[, k], breaks=seq(-5, 5, 0.1),
     main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
     xlab="theta", yaxt='n', ylab="")
abline(v=theta_true[k], col="white", lwd=2)
abline(v=theta_true[k], col="black", lwd=1)

# Hierarchical
fit <- stan(file='hier_cauchy_cp.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

samples = extract(fit)

par(mfrow=c(1, 1))

hist(abs(rnorm(4000, 0, 10)), breaks=seq(0, 50, 0.5),
     main="", col=c_light, border=c_light_highlight,
     xlab="tau", yaxt='n', ylab="", ylim=c(0, 800))
hist(samples$tau, breaks=seq(0, 50, 0.5),
     main="", col=c_dark, border=c_dark_highlight, add=T)

abline(v=0.1, col="white", lwd=4)
abline(v=0.1, col="black", lwd=2)

abline(v=20, col="white", lwd=4)
abline(v=20, col="black", lwd=2)

par(mfrow=c(3, 3))

for (k in 1:9) {
  hist(samples$theta[, k], breaks=seq(-40, 40, 0.25),
       main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
       xlab="theta", yaxt='n', ylab="")
  abline(v=theta_true[k], col="white", lwd=2)
  abline(v=theta_true[k], col="black", lwd=1)
}

par(mfrow=c(1, 1))

k <- 1
hist(samples$theta[, k], breaks=seq(-5, 5, 0.1),
     main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
     xlab="theta", yaxt='n', ylab="")
abline(v=theta_true[k], col="white", lwd=2)
abline(v=theta_true[k], col="black", lwd=1)

############################################################
#
# One Horse Town 
#
############################################################

# Before fitting the full horseshoe population model let's fix the 
# population scale.  A shoe for geldings, if you will.

data$tau <- 0.1

fit <- stan(file='geldingshoe_cp.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

# Eeeek, divergences and some indications of metastability
# in the Markov chains.  Let's go hunting.

partition <- util$partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mfrow=c(3, 3))

for (k in 1:9) {
  name_x <- paste("theta[", k, "]", sep='')
  name_y <- paste("lambda[", k, "]", sep='')
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("k = ", k),
       xlab=name_x, xlim=c(theta_true[k] - 5, theta_true[k] + 5), 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-4, 8))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
  abline(v=alpha_true[k], col="white", lwd=2)
  abline(v=alpha_true[k], col="black", lwd=1)
  abline(h=0, col="white", lwd=2)
  abline(h=0, col="black", lwd=1)
}

# Lots of funnels between theta and lambda as expected from the 
# hierarchical coupling!  Because each lambda is informed by only 
# a single context there's not much to fight against the funnel 
# geometry in the centered population model.

# Note that the funnels in contexts 7, 8, and 9 push a little bit higher.
# Interestingly these are the contexts with the largest values of theta
# in their true data generating processes.  

# This is a fundamental property of the horseshoe population model.
# When lambda is around one the posterior distribution for the
# corresponding theta is strongly influenced by the inner core,
# favoring a non-centered parameterization even when the individual 
# likelihood functions are reasonably well-informed.  On the other 
# hand when lambda is much larger than one the horseshoe model
# decouples from the core prior model and the posterior distribution 
# for the corresponding theta is dominated by the individual likelihood 
# function.

# Perhaps we can employ a mixed parameterization, centering just those 
# contexts where theta is inferred to be much larger than the others.

data$cp_idx <- c(7, 8, 9)
data$K_cp <- length(data$cp_idx)

data$ncp_idx <- setdiff(1:data$K, data$cp_idx)
data$K_ncp <- length(data$ncp_idx)

fit <- stan(file='geldingshoe_mixed.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

# Didn't seem to have done much.  Indeed now the inverted funnels 
# in the other contexts induce all of the divergent transitions.

partition <- util$partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mfrow=c(3, 3))

for (k in 1:9) {
  if (k %in% data$cp_idx)
    name_x <- paste("theta_cp[", which(data$cp_idx == k), "]", sep='')
  else
    name_x <- paste("theta_ncp[", which(data$ncp_idx == k), "]", sep='')
  
  name_y <- paste("lambda[", k, "]", sep='')
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("k = ", k),
       xlab=name_x, xlim=c(theta_true[k] - 5, theta_true[k] + 5), 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-7, 8))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
}

# Unfortunately the Cauchy prior on lambda which allows for the
# horseshoe model to decouple from the core population also 
# amplifies the funnel geometry in both parameterizations.  The
# heavy tail not only pulls the funnel deeper but also increases 
# how quickly the curvature changes.  We've taken funnels and 
# pulled them into more pathological geometries!

# With both the centered and non-centered parameterizations of the 
# strongly regularized contexts exhibiting funnel geometries let's 
# pull out the heavy artillery and consider a partially-centered
# parameterization that interpolates between the two.

data$w <- c(0.5, 0.5, 0.5, 
            0.5, 0.5, 0.5,
            1.0, 1.0, 1.0)

fit <- stan(file='geldingshoe_partial.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

# Not much help.  

# Looks like these partially-centered parameterizations are too
# centered.

partition <- util$partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mfrow=c(3, 3))

for (k in 1:9) {
  name_x <- paste("theta_tilde[", k, "]", sep='')
  name_y <- paste("lambda[", k, "]", sep='')
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("k = ", k),
       xlab=name_x, xlim=c(theta_true[k] - 5, theta_true[k] + 5), 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-7, 8))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
}

# Let's try a more non-centered parameterization.

data$w <- c(0.25, 0.25, 0.25, 
            0.25, 0.25, 0.25,
            1.0, 1.0, 1.0)

fit <- stan(file='geldingshoe_partial.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

partition <- util$partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mfrow=c(3, 3))

for (k in 1:9) {
  name_x <- paste("theta_tilde[", k, "]", sep='')
  name_y <- paste("lambda[", k, "]", sep='')
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("k = ", k),
       xlab=name_x, xlim=c(theta_true[k] - 5, theta_true[k] + 5), 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-7, 8))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
}

# Even more non-centered!

data$w <- c(0.2, 0.2, 0.2, 
            0.2, 0.2, 0.2,
            1.0, 1.0, 1.0)

fit <- stan(file='geldingshoe_partial.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

# We've found a nice balance between regular and inverted 
# funnels, but that just means that they're both pathological 
# now.  Cauchy density functions are the worst.

partition <- util$partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mfrow=c(3, 3))

for (k in 1:9) {
  name_x <- paste("theta_tilde[", k, "]", sep='')
  name_y <- paste("lambda[", k, "]", sep='')
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("k = ", k),
       xlab=name_x, xlim=c(theta_true[k] - 5, theta_true[k] + 5), 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-7, 8))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
}

# Maybe we can brute force it by pushing
# Stan's adaptation to be much less aggressive than default.

fit <- stan(file='geldingshoe_partial.stan', data=data, 
            seed=4938483, refresh=1000, control=list(adapt_delta=0.99))

util$check_all_diagnostics(fit)

# The marginal posteriors for each context parameter cover the
# true values, both for the irrelevant values and the relevant
# values

samples = extract(fit)

par(mfrow=c(9, 2))

for (k in 1:9) {
  hist(samples$theta[, k], breaks=seq(theta_true[k] - 10, theta_true[k] + 10, 0.25),
       main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
       xlab="theta", yaxt='n', ylab="")
  abline(v=theta_true[k], col="white", lwd=2)
  abline(v=theta_true[k], col="black", lwd=1)
  
  hist(log(samples$lambda[, k]), breaks=seq(-10, 10, 0.25),
       main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
       xlab="log(lambda)", yaxt='n', ylab="")
  abline(v=0, col="white", lwd=2)
  abline(v=0, col="black", lwd=1)
}


############################################################
#
# High Horse 
#
############################################################

# Let the scale of the inner core population be free.

# Using what we learned before let's use a mixed parameterization again.

data$cp_idx <- c(7, 8, 9)
data$K_cp <- length(data$cp_idx)

data$ncp_idx <- setdiff(1:data$K, data$cp_idx)
data$K_ncp <- length(data$ncp_idx)

fit <- stan(file='horseshoe_mixed.stan', data=data, 
             seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

# Didn't seem to have done much.  Indeed now the inverted funnels 
# in the other contexts induce all of the divergent transitions.

partition <- util$partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mfrow=c(3, 3))

for (k in 1:9) {
  if (k %in% data$cp_idx)
    name_x <- paste("theta_cp[", which(data$cp_idx == k), "]", sep='')
  else
    name_x <- paste("theta_ncp[", which(data$ncp_idx == k), "]", sep='')

  name_y <- paste("lambda[", k, "]", sep='')
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("k = ", k),
       xlab=name_x, 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-9, 9))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
  abline(h=0, col="white", lwd=2)
  abline(h=0, col="black", lwd=1)
}

# Maybe we can brute force it.  Because the centered funnels were 
# a bit less pathological let's center everything but also push 
# Stan's adaptation to be much less aggressive than default.

fit <- stan(file='horseshoe_mixed.stan', data=data, 
             seed=4938483, refresh=1000, control=list(adapt_delta=0.99))

util$check_all_diagnostics(fit)

# This does reduce the divergent transitions, although sadly 
# it does not eliminate them.  The more accurate fit, however,
# does better show just how narrow those Cauchy funnels are.

partition <- util$partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mfrow=c(3, 3))

for (k in 1:9) {
  if (k %in% data$cp_idx)
    name_x <- paste("theta_cp[", which(data$cp_idx == k), "]", sep='')
  else
    name_x <- paste("theta_ncp[", which(data$ncp_idx == k), "]", sep='')
  
  name_y <- paste("lambda[", k, "]", sep='')
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("k = ", k),
       xlab=name_x, 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-8, 8))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
  abline(h=0, col="white", lwd=2)
  abline(h=0, col="black", lwd=1)
}

# Can we save the horseshoe by partially centering?

data$w <- c(0.5, 0.5, 0.5, 
            0.5, 0.5, 0.5,
            1.0, 0.5, 1.0)

fit <- stan(file='horseshoe_partial.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

# Looks like these partially-centered parameterizations are too
# centered.  Note also that some of the contexts show signs of 
# funnels bothon top and on bottom at the same time!

partition <- util$partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mfrow=c(3, 3))

for (k in 1:9) {
  name_x <- paste("theta_tilde[", k, "]", sep='')
  name_y <- paste("lambda[", k, "]", sep='')
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("Context", k),
       xlab=name_x, 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-8, 8))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
  abline(h=0, col="white", lwd=2)
  abline(h=0, col="black", lwd=1)
}

# Too centered.  Let's go more non-centered.

data$w <- c(0.25, 0.25, 0.25, 
            0.25, 0.25, 0.25,
            1.00, 0.25, 1.00)

fit <- stan(file='horseshoe_partial.stan', data=data, 
            seed=4938483, refresh=1000)

util$check_all_diagnostics(fit)

# Shine bright like a diamond!

partition <- util$partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mfrow=c(3, 3))

for (k in 1:9) {
  name_x <- paste("theta_tilde[", k, "]", sep='')
  name_y <- paste("lambda[", k, "]", sep='')
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("Context", k),
       xlab=name_x, 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-7, 4))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
  abline(h=0, col="white", lwd=2)
  abline(h=0, col="black", lwd=1)
}

# Finally let's try relaxing the adaptation for this partial
# centering.

fit <- stan(file='horseshoe_partial.stan', data=data, 
            seed=4938483, refresh=1000, 
            control=list(adapt_delta=0.999, max_treedepth=15))

util$check_all_diagnostics(fit)

# We're exhausted.  Let's ignore that lingering divergences and
# take a look at the (mildly biased) posterior inferences.

samples = extract(fit)

# The marginal posteriors for each context parameter cover the
# true values, both for the irrelevant values and the relevant
# values.

par(mfrow=c(9, 2))

for (k in 1:9) {
  hist(samples$theta[, k], breaks=seq(theta_true[k] - 10, theta_true[k] + 10, 0.25),
       main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
       xlab="alpha", yaxt='n', ylab="")
  abline(v=theta_true[k], col="white", lwd=2)
  abline(v=theta_true[k], col="black", lwd=1)
  
  hist(log(samples$lambda[, k]), breaks=seq(-10, 10, 0.25),
       main=paste("k = ", k), col=c_dark, border=c_dark_highlight,
       xlab="log(lambda)", yaxt='n', ylab="")
  abline(v=0, col="white", lwd=2)
  abline(v=0, col="black", lwd=1)
}

par(mfrow=c(1, 1))

k <- 1
hist(rnorm(4000, theta_true[k], 1), breaks=seq(-5, 5, 0.1),
     main=paste("k = ", k), col=c_light, border=c_light_highlight,
     xlab="alpha", yaxt='n', ylab="", ylim=c(0, 1250))
hist(samples$theta[, k], breaks=seq(-5, 5, 0.1),
     col=c_dark, border=c_dark_highlight, add=T)
abline(v=theta_true[k], col="white", lwd=2)
abline(v=theta_true[k], col="black", lwd=1)

# At least we have some comfort that we recover the true
# relevance scale reasonably well, especially given that we
# have only nine contexts across which we can partially pool.

par(mfrow=c(1, 1))

hist(abs(rnorm(4000, 0, 0.1)), breaks=seq(0, 0.5, 0.005),
     main="", col=c_light, border=c_light_highlight,
     xlab="tau", yaxt='n', ylab="", ylim=c(0, 200))
hist(samples$tau, breaks=seq(0, 0.5, 0.005),
     main="", col=c_dark, border=c_dark_highlight, add=T)
abline(v=1, col="white", lwd=2)
abline(v=0.1, col="white", lwd=2)
abline(v=0.1, col="black", lwd=1)

# In some sense the horseshoe population model works as advertised.
# The individual lambda give the model freedom to decouple large 
# parameters from the prior regularization while pooling the rest 
# towards the irrelevant population core.

# That decoupling, however, comes at the cost of Cauchy funnels 
# that frustrate sampler exploration and are nearly impossible 
# to reparameterize away.  Even knowing which parameters were 
# large and which were small we still had to exploit a carefully 
# tuned partially centered parameterization and push down to a 
# small step size, and expensive numerical trajetories, in order 
# explore the full horseshoe geometry.

