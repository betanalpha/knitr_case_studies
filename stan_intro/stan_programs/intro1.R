############################################################
# Initial Setup
############################################################

library(rstan)
rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores()) # Parallelize chains

util <- new.env()
source('stan_utility.R', local=util)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

############################################################
# Normal Model
############################################################

# Our first Stan program needs the number of ys in the form of a list
N <- 3
data <- list("N" = N)

# Compile Stan program and fit with Hamiltonian Monte Carlo
fit <- stan(file='normal1.stan', data=data, seed=4938483)

# Check diagnostics one by one
util$check_n_eff(fit)
util$check_rhat(fit)
util$check_div(fit)
util$check_treedepth(fit)
util$check_energy(fit)

# Or all at once
util$check_all_diagnostics(fit)

# Plot marginal distributions
params <- extract(fit)

par(mfrow=c(2, 2))
   
hist(params$y[,1], breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="y[1]", yaxt='n', ylab="")
     
hist(params$y[,2], breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="y[2]", yaxt='n', ylab="")

hist(params$y[,3], breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="y[3]", yaxt='n', ylab="")

hist(params$theta, breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="theta", yaxt='n', ylab="")
     
# Let's repeat using the vectorized Stan program
fit <- stan(file='normal2.stan', data=data, seed=4938483)

util$check_all_diagnostics(fit)

params <- extract(fit)

par(mfrow=c(2, 2))
   
hist(params$y[,1], breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="y[1]", yaxt='n', ylab="")

hist(params$y[,2], breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="y[2]", yaxt='n', ylab="")

hist(params$y[,3], breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="y[3]", yaxt='n', ylab="")

hist(params$theta, breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="theta", yaxt='n', ylab="")
     
# And finally the vectorized Stan program using sampling
# statements instead of incrementing target directly
fit <- stan(file='normal3.stan', data=data, seed=4938483)

util$check_all_diagnostics(fit)

params <- extract(fit)

par(mfrow=c(2, 2))
  
hist(params$y[,1], breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="y[1]", yaxt='n', ylab="")

hist(params$y[,2], breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="y[2]", yaxt='n', ylab="")

hist(params$y[,3], breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="y[3]", yaxt='n', ylab="")

hist(params$theta, breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="theta", yaxt='n', ylab="")

############################################################
# Generated Quantities
############################################################

# What if we wanted to compute the expectation of exp(theta)
# in addition to identity(theta)?  We can accomplish this by
# either adding exp(theta) to the "transformed parameters"
# block or the "generated quantities" block.  Becuase we
# don't need exp(theta) in the model block the latter will
# be the best choice.

fit <- stan(file='normal_gq.stan', data=data, seed=4938483)

util$check_all_diagnostics(fit)

params <- extract(fit)

par(mfrow=c(2, 3))

hist(params$y[,1], breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="y[1]", yaxt='n', ylab="")

hist(params$y[,2], breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="y[2]", yaxt='n', ylab="")

hist(params$y[,3], breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="y[3]", yaxt='n', ylab="")

hist(params$theta, breaks=seq(-7, 7, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(-7, 7), xlab="theta", yaxt='n', ylab="")

hist(params$exp_theta, breaks=seq(0, 60, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, main="",
     xlim=c(0, 60), xlab="exp_theta", yaxt='n', ylab="")

############################################################
# Conditional Normal Model
############################################################
     
# Instead of fitting the ys jointly with theta let's
# condition the joint probability density on some
# specific values of the ys
  
N <- 3
y <- c(-1, -0.25, 0.3)
data <- list("N" = N, "y" = y)
     
# Compile Stan program and fit with Hamiltonian Monte Carlo
cond_fit <- stan(file='conditional_normal.stan', data=data, seed=4938483)

# Check diagnostics
util$check_all_diagnostics(cond_fit)

# Plot conditional distribution
cond_params <- extract(cond_fit)

par(mfrow=c(1, 2))

hist(params$theta, breaks=seq(-4, 4, 0.25), prob=T,
     col=c_light, border=c_light_highlight, main="",
     xlim=c(-4, 4), xlab="theta", yaxt='n', ylab="", ylim=c(0, 0.75))
hist(cond_params$theta, breaks=seq(-4, 4, 0.25), prob=T,
     col=c_dark, border=c_dark_highlight, add=T)

hist(params$exp_theta, breaks=seq(0, 60, 0.6), prob=T,
     col=c_light, border=c_light_highlight, main="",
     xlim=c(0, 60), xlab="exp_theta", yaxt='n', ylab="")
hist(cond_params$exp_theta, breaks=seq(0, 60, 0.5), prob=T,
     col=c_dark, border=c_dark_highlight, add=T)
