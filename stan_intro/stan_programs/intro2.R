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

c_light_trans <- c("#DCBCBC80")
c_light_highlight_trans <- c("#C7999980")
c_mid_trans <- c("#B97C7C80")
c_mid_highlight_trans <- c("#A2505080")
c_dark_trans <- c("#8F272780")
c_dark_highlight_trans <- c("#7C000080")

############################################################
# Gamma Model
############################################################
    
# What happens when we try to fit with an unconstrained parameter?
fit <- stan(file='gamma1.stan', seed=4938483)

# What's with all of the divergences?!
util$check_all_diagnostics(fit)

params <- extract(fit)

theta <- seq(0, 8, 0.001)
plot(theta, dgamma(theta, 1.25, 1.25), type="l", col=c_dark_highlight, 
    lwd=2, xlab="theta", ylab="Probability Density", yaxt='n')
hist(params$theta, breaks=seq(0, 8, 0.1),
     col=c_mid_trans, border=c_mid_highlight_trans, probability=T, add=T)

# We forget to constraint theta to be positive in the Stan
# program!  Let's add the appropriate constraint and try again
fit <- stan(file='gamma2.stan', seed=4938483)

util$check_all_diagnostics(fit)

params <- extract(fit)

theta <- seq(0, 8, 0.001)
plot(theta, dgamma(theta, 1.25, 1.25), type="l", col=c_dark_highlight, 
     lwd=2, xlab="theta", ylab="Probability Density", yaxt='n')
hist(params$theta[params$theta < 8], breaks=seq(0, 8, 0.1),
     col=c_mid_trans, border=c_mid_highlight_trans, probability=T, add=T)

############################################################
# Truncated Model
############################################################

# What happens when we don't properly normalize a truncated model?
fit1 <- stan(file='truncation1.stan', seed=4938483)
util$check_all_diagnostics(fit1)

fit2 <- stan(file='truncation2.stan', seed=4938483, control=list(adapt_delta=0.85))
util$check_all_diagnostics(fit2)

hist(extract(fit1)$x, breaks=seq(-5, 5, 0.25), 
     main='', xlab="x", ylab="Probability Density", yaxt='n', 
     col=c_dark_trans, border=c_dark_highlight_trans, probability=T)
hist(extract(fit2)$x, breaks=seq(-5, 5, 0.25), 
     col=c_light_trans, border=c_light_highlight_trans, probability=T, add=T)

# The normalization, however, doesn't matter if the trunctation is known
fit3 <- stan(file='truncation3.stan', seed=4938483, control=list(adapt_delta=0.85))
util$check_all_diagnostics(fit3)

fit4 <- stan(file='truncation4.stan', seed=4938483, control=list(adapt_delta=0.85))
util$check_all_diagnostics(fit4)

hist(extract(fit3)$x, breaks=seq(-5, 5, 0.25), 
     main='', xlab="x", ylab="Probability Density", yaxt='n', 
     col=c_dark_trans, border=c_dark_highlight_trans, probability=T)
hist(extract(fit4)$x, breaks=seq(-5, 5, 0.25), 
     col=c_light_trans, border=c_light_highlight_trans, probability=T, add=T)
