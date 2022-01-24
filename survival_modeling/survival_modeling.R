##################################################
# Colors
##################################################

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_trans <- c("#DCBCBC80")
c_dark_trans <- c("#8F272780")
c_green_trans <- c("#00FF0080")

library(colormap)
nom_colors <- c("#DCBCBC", "#C79999", "#B97C7C", "#A25050", "#8F2727", "#7C0000")

par(family="CMU Serif", las=1, bty="l", cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 5))

##################################################
# Set up Stan
##################################################

library(rstan)
rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores()) # Parallelize chains
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

util <- new.env()
source('stan_utility.R', local=util)

############################################################
# Mathematical Setup
############################################################

# Consider an event that occurs after t = 0 but before t = infty.
# Let pi(t) be the density representing the distribution of times
# for the occurance of an event that can begin as early as t = 0.
# Assuming that nothing about the system is observed other than the
# occurance of the event itself.

# The cumulative distribution function quantifies the probability
# that an event will occur before some time t,
#   Pi(t) = \int_0^t dt' pi(t') = pi(occurance before t)
# Likewise the complementary cumulative distribution function
# quantifies the probability that na event has not occured before
# some time t,
#   S(t) = 1 - Pi(t) = \int_t^infty dt' pi(t') = pi(no occurance before t)

# Modeling pi(t), however, isn't always straightforward.  Often more
# straightforward to model the rate of occurance after _waiting_ for 
# some time.  To that end we can introduce the probability density
# for an event to occur immediately only after waiting until t,
#  pi_wait(t) = pi(t | no occurance before t).

# The probability density function for this monitored event occurance is 
# given by
#   pi(t) = pi(t | no occurance before t) * pi(no occurance before t)
#   pi(t) = pi_wait(t) * int_t^infty dt' pi(t')
#         = pi_wait(t) * S(t)
# or
#  pi_wait(t) = pi(t) / S(t).

# This provides an explicit result for pi_wait(t) given pi(t).  If we model 
# pi_wait(t), however, this provides only an implicit form.  To give an 
# explicit form we take advantage of the fact that
#  dS/dt(t) = d/dt (1 - \int_0^t dt' pi(t'))
#            = - d/dt \int_0^t dt' pi(t')
#            = - pi(t')
# In other words our implicit definition of pi(t) defines a differential 
# equation for the survival function,
#   pi(t) / S(t) = pi_wait(t)
#   - dS/dt(t) / S(t) = pi_wait(t)
#   1/S(t) dS/dt(t) = - pi_wait(t)
#   dlog S/dt (t) = - pi_wait(t).
# Using the initial condition S(0) = 1, or log S(0) = 0, we can integrate
# this differential equation to give
#   log S(t) = - \int_0^{t} dt' pi_wait(t')
#   S(t) = exp( - \int_0^{t} dt' pi_wait(t') )

# Be careful with the multiple integrals here:
#  1 - \int_0^t dt' pi(t') = S(t) = exp( - \int_0^{t} dt' pi_wait(t') )

# Substituting this back into the definition for pi(t) gives
#   pi(t) = pi_wait(t) * exp( - \int_0^{t} dt' pi_wait(t') )
# or equivalently
#   log pi(t) = log pi_wait(t) - \int_0^t dt' pi_wait(t')
# Given pi_wait(t) we need to be able to integrate to give pi.
# In some cases this can be done analytically and in others we have to
# rely on numerical integration to evaluate pi(t) at each time of
# interest.

############################################################
# The Hazard Function
############################################################

# Typical survival models start with an assumed form for pi_wait(t),
# which is called the _instantateous hazard function_, or _hazard function_ 
# for short.  Because pi_wait(t) is a conditional density function
# is cannot be negative.  At the same time because we're conditioning
# on a given time pi_wait(t) is not a well-defined denstiy function over 
# all of t at once.  In particular because S(infty) = 0 we must have
#   0 = S(\infty) = exp( - \int_0^{infty} dt' pi_wait(t') )
# or
#   \int_0^{infty} dt' pi_wait(t') = infty.
# In other words pi_wait(t) isn't normalizeable over all of [0, infty).

# The log survival function is given by integrating all of the instantatenous 
# hazards together into a cumulative hazard function, 
#  log S(t) = - int_0^t dt' pi_wait(t') = Lambda(t).
# The more hazards that have been accumulated up to a given time the smaller
# the surival function and the more likely the event is to have occurred
# by that time.

# The hazard function quantifies instantaenous changes to the survival 
# function over time.  If the hazard function is zero for some interval, 
# for example, then the survival function is constant across that interval.
# Positive hazards imply decreasing survival, and the larger the hazard
# the faster the survival function decreases.  Because the hazard 
# function is positive, however, the accumulated hazard can never 
# decrease.

############################################################
# Simulating Survival Events
############################################################

# Once we've constructed pi(t) simulation of uncensored and 
# censored events is straightforward, at least in theory.

# For uncensored events we just sample from pi(t).  If we 
# don't have access to an appropriate random number generator
# then we may be able to use the inverse CDF method,
#   p ~ U(0, 1)
#   t = S^{-1}(p),
# so long as we can invert the survival function.

############################################################
# Exponential Model
############################################################

# Constant hazard, $pi_wait(t) = lambda$.
# Well-defined because conditioned on no previous occurance?

# Cumulative hazard = Lambda(t) = int_0^t lambda = lambda * t.

# Survival function = S(t) = exp(-Lambda(t)) = exp(-lambda * t)

# Final Density = pi_wait(t) * S(t) = lambda * exp(-lambda * t) = exponential(t | lambda).

# Implement with survival model and exponential model directly.

# For simulations:
# p = exp(-lambda * t)
# log p = - lambda * t
# t = - log p / lambda.

# Visualize True Model Configuration
gamma <- 0.1

M <- 100
xs <- seq(0, 26, 26 / (M - 1))
ys <- dexp(xs, gamma)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylab="Probability Density Function", ylim=c(0, 0.2))

ts <- rexp(20, 0.2)
cat(sprintf("%.3f,", sort(ts), "\n"))
print("")


M <- 100
xs <- seq(0, 2, 2 / (M - 1))

par(mar = c(5, 4, 3, 2))
par(mfrow=c(2, 2))

ys <- rep(gamma, M)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 6), ylab="Hazard Function")

ys <- gamma * xs
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 6), ylab="Cumulative Hazard Function")

ys <- exp(-gamma * xs)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 1), ylab="Survival Function")

ys <- gamma * exp(-gamma * xs)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylab="Probability Density Function", yaxt='n')

# Simulate data
simu_data <- list("N" = 100, "gamma" = gamma)

simu <- stan(file='simu_const_haz.stan', data=simu_data,
             iter=1, chains=1, seed=4838282, 
             algorithm="Fixed_param")
obs_times <- array(extract(simu)$obs_times[1,])

data <- list("N" = simu_data$N, "obs_times" = obs_times)

# Visualize Data
plot_eccdf <- function(obs_times) {
  N <- length(obs_times)
  
  ordered_data <- sort(obs_times)
  
  xs <- c(ordered_data[1] - 0.5,
          rep(ordered_data, each=2),
          ordered_data[N] + 0.5)

  eccdf <- rep(N:0, each=2) / N

  lines(xs, eccdf, lwd="3", col="white")
  lines(xs, eccdf, lwd="2", col="black")
}

par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))

M <- 100
xs <- seq(0, 2, 2 / (M - 1))
ys <- exp(-gamma * xs)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 1), ylab="Survival")

plot_eccdf(data$obs_times)

# Plot inferred survival behavior
data$N_display <- 100
data$display_times <- seq(0, 2, 2 / (data$N_display - 1))

# Fit with survival model
fit <- stan(file='fit_const_haz.stan', data=data, 
            seed=4938483, refresh=0)

# False positives due to constant survival at t = 0
util$check_all_diagnostics(fit)

surv_samples = extract(fit)

plot_marginal <- function(name, posterior_samples, truth,
                          display_xlims, title="") {
  posterior_values <- posterior_samples[name][[1]]
  bin_lims <- range(posterior_values)
  delta <- diff(range(posterior_values)) / 50
  breaks <- seq(bin_lims[1], bin_lims[2]+ delta, delta)
        
  hist(posterior_values, breaks=breaks,
       main=title, xlab=name, xlim=display_xlims,
       ylab="", yaxt='n',
       col=c_dark, border=c_dark_highlight)
        
  abline(v=truth, col="white", lty=1, lwd=3)
  abline(v=truth, col="black", lty=1, lwd=2)
}

par(mar = c(5, 2, 3, 2))
par(mfrow=c(1, 1))

plot_marginal("gamma", surv_samples, gamma, c(0, 5))

plot_marginal_quantiles <- function(grid, f_samples, xlab, xlim, ylab, ylim) {
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:length(grid),
                 function(n) quantile(f_samples[,n], probs=probs))
        
  plot(1, type="n", main="",
       xlab=xlab, xlim=xlim, ylab=ylab, ylim=ylim)
  polygon(c(grid, rev(grid)), c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(grid, rev(grid)), c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(grid, rev(grid)), c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(grid, rev(grid)), c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
        
  lines(grid, cred[5,], col=c_dark, lwd=2)
}

par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))
plot_marginal_quantiles(data$display_times, surv_samples$display_survival,
                        xlab="Time", xlim=c(0, 2),
                        ylab="Survival", ylim=c(0, 1))
plot_eccdf(data$obs_times)

plot_hist_retro <- function(data, pred_samples, name, min, max, delta,
                            xlim_override=NA) {
        
  breaks <- seq(min, max, delta)
  B <- length(breaks) - 1
        
  idx <- rep(1:B, each=2)
  xs <- sapply(1:length(idx),
               function(b) if(b %% 2 == 1) breaks[idx[b]] else breaks[idx[b] + 1])
        
  obs <- hist(data, breaks=breaks, plot=FALSE)$counts
  pad_obs <- do.call(cbind, lapply(idx, function(n) obs[n]))
        
  post_pred <- sapply(1:4000,
                      function(n) hist(pred_samples[n,], breaks=breaks, plot=FALSE)$counts)
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:B, function(b) quantile(post_pred[b,], probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9, n]))
        
  xlim <- c(min, max)
  if (!is.na(xlim_override))
    xlim <- xlim_override
        
  plot(1, type="n", main="Posterior Retrodictive Check",
       xlim=xlim, xlab=name,
       ylim=c(0, max(c(obs, cred[9,]))), ylab="")
        
  polygon(c(xs, rev(xs)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  for (b in 1:B)
    lines(xs[(2 * b - 1):(2 * b)], pad_cred[5,(2 * b - 1):(2 * b)], 
          col=c_dark, lwd=2)
        
  lines(xs, pad_obs, col="white", lty=1, lw=2.5)
  lines(xs, pad_obs, col="black", lty=1, lw=2)
  
  post_pred <- sapply(1:4000,
                      function(n) hist(pred_samples[n,], breaks=breaks, plot=FALSE)$counts 
                                  - obs)
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:B, function(b) quantile(post_pred[b,], probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9, n]))
  
  xlim <- c(min, max)
  if (!is.na(xlim_override))
          xlim <- xlim_override
  
  plot(1, type="n", main="Residuals",
       xlim=xlim, xlab=name,
       ylim=c(min(cred[1,]), max(cred[9,])), ylab="Posterior Predictive - Observations")
  
  abline(h=0, col="#DDDDDD", lty=2, lwd=2)
  
  polygon(c(xs, rev(xs)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  for (b in 1:B)
          lines(xs[(2 * b - 1):(2 * b)], pad_cred[5,(2 * b - 1):(2 * b)], 
                col=c_dark, lwd=2)
}

par(mfrow=c(1, 2))
plot_hist_retro(data$obs_times, surv_samples$pred_obs_times, 
                "Event Times", 0, 5, 0.25, c(0, 2))

# Fit with exponential model
fit <- stan(file='fit_exp.stan', data=data, 
            seed=4938483, refresh=0)

util$check_all_diagnostics(fit)

exp_samples = extract(fit)


par(mar = c(5, 2, 3, 2))
par(mfrow=c(1, 2))

plot_marginal("gamma", surv_samples, gamma, c(0, 5), 
              "Constant Hazard Survival Model")
plot_marginal("gamma", exp_samples, gamma, c(0, 5),
              "Exponential Model")

############################################################
# Weibull Model
############################################################

# lambda(t) = alpha * (1/sigma) * (t / sigma)^{alpha - 1}

# Lambda(t) = (t / sigma)^alpha

# f(t) = lambda(t) * exp(-Lambda(t))
#      = alpha * (1/sigma) * (t / sigma)^{alpha - 1} * exp( - (t / sigma)^alpha )

# Implement with survival model and Weibull model directly.
        
# For simulations:
# p = exp(- (t / sigma)^alpha)
# log p = - (t / sigma)^alpha
# t / sigma = (- log p)^{1 / alpha}
# t = sigma * (- log p)^{1 / alpha}

# Visualize True Model Configuration
alpha <- 1.25
sigma <- 2

M <- 100
xs <- seq(0, 10, 10 / (M - 1))

par(mar = c(5, 4, 3, 2))
par(mfrow=c(2, 2))

ys <- (alpha / sigma) * (xs / sigma)**(alpha - 1)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 8), ylab="Hazard Function")

ys <- (xs / sigma)**alpha
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 8), ylab="Cumulative Hazard Function")

ys <- exp(-(xs / sigma)**alpha)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 1), ylab="Survival Function")

ys <- (alpha / sigma) * (xs / sigma)**(alpha - 1) * exp(-(xs / sigma)**alpha)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylab="Probability Density Function", ylim=c(0, 0.4))
        
# Simulate data
simu_data <- list("N" = 100, "alpha" = alpha, "sigma" = sigma)

simu <- stan(file='simu_pow_haz.stan', data=simu_data,
             iter=1, chains=1, seed=4838282, 
             algorithm="Fixed_param")
obs_times <- array(extract(simu)$obs_times[1,])

data <- list("N" = simu_data$N, "obs_times" = obs_times)

# Visualize Data
par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))

M <- 100
xs <- seq(0, 8, 8 / (M - 1))
ys <-  exp(-(xs / sigma)**alpha)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 1), ylab="Survival")

plot_eccdf(data$obs_times)

# Plot inferred survival behavior
data$N_display <- 100
data$display_times <- seq(0, 8, 8 / (data$N_display - 1))

# Fit with survival model
fit <- stan(file='fit_pow_haz.stan', data=data, 
            seed=4938483, refresh=0)

# False positives due to constant survival at t = 0
util$check_all_diagnostics(fit)

surv_samples = extract(fit)

par(mar = c(5, 2, 3, 2))
par(mfrow=c(1, 2))

plot_marginal("alpha", surv_samples, alpha, c(0, 5))
plot_marginal("sigma", surv_samples, sigma, c(0, 5))

par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))
plot_marginal_quantiles(data$display_times, surv_samples$display_survival,
                        xlab="Time", xlim=c(0, 8),
                        ylab="Survival", ylim=c(0, 1))
plot_eccdf(data$obs_times)

par(mfrow=c(1, 2))
plot_hist_retro(data$obs_times, surv_samples$pred_obs_times, 
                "Event Times", 0, 10, 0.5, c(0, 8))

# Fit with Weibull model
fit <- stan(file='fit_weibull.stan', data=data, 
            seed=4938483, refresh=0)

util$check_all_diagnostics(fit)

weibull_samples = extract(fit)

par(mar = c(5, 2, 3, 2))
par(mfrow=c(2, 2))

plot_marginal("alpha", surv_samples, alpha, c(0, 5), 
              "Power Hazard Survival Model")
plot_marginal("alpha", weibull_samples, alpha, c(0, 5),
              "Weibull Model")

plot_marginal("sigma", surv_samples, sigma, c(0, 5),
              "Power Hazard Survival Model")
plot_marginal("sigma", weibull_samples, sigma, c(0, 5),
              "Weibull Model")

############################################################
# Censoring
############################################################

# In practice, however, we can't monitor a system forever.  Most survival 
# experiments continue only to some maximum time; if the event 
# has not yet occured by this time then it becomes censored.  
# Fortunately the censored model is given exactly by the
# surival function,
#  pi(t > t_max) = \int_{t_max}^infty dt' pi(t') = S(t_max).
# These events are referred to as _right-censored_ events.

# If we begin monitoring only at some time t_min > 0 then any events
# that occur before t_min are similarly _left-censored_,
#  pi(t < t_min) = \int_0^{t_min} dt' pi(t') = 1 - S(t_min)

# An event that is both left-censored and right-censored is 
# referred to as an _interval-censored_ event and is modeled with
# the probability
#   pi(t_l < t < t_u) = \int_{t_l}^{t_u} dt' pi(t')
#                          =   \int_{t_l}^infty dt' pi(t')
#                            - \int_{t_u}^infty dt' pi(t')
#                          = S(t_l) - S(t_u)

# For censored events we categorize all of the possible observation
# intervals, sample from the corresponding categorical to determine
# into which interval the simulated event falls, and then sample any 
# necessary structure within that simulated interval.

# Or just simulate and then bin based on values

# Sample from events with times (t1, t2)
#   p ~ U(S(t2), S(t1))
#   t = S^{-1}(p).

# For example if we have right-censoring then we need to compute
# pi(t <= t_max) = 1 - S(t_max) and pi(t > t_max) = S(t_max) 
# and sample from the corresponding Bernoulli.  If we sample t > t_max 
# then we're done, but if we sample t <= t_max then we need to sample 
# from the truncated waiting time distribution.  When we can invert
# the survival function this is straightforward,
#   p ~ U(S(t_infty), S(t_max))
#   t = S^{-1}(p).
# equivalent to
#  p ~ U(0, S(t_max))
#  t = S^{-1}(p)

# Example with interval censoring.

# Example with all of the censorings.

# If we can't work out the survival function in closed form, let 
# alone invert it, then we have to work out t = S^{-1}(p) numerically,
# for example with a root finding method.

# Visualize True Model Configuration
gamma <- 0.5

M <- 100
xs <- seq(0, 10, 10 / (M - 1))

par(mar = c(5, 4, 3, 2))
par(mfrow=c(2, 2))

ys <- rep(gamma, M)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 6), ylab="Hazard Function")

ys <- gamma * xs
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 6), ylab="Cumulative Hazard Function")

ys <- exp(-gamma * xs)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 1), ylab="Survival Function")

ys <- gamma * exp(-gamma * xs)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylab="Probability Density Function", yaxt='n')

# Visualize censoring intervals
t1 <- 2
t2 <- 6
t3 <- 8

par(mar = c(5, 2, 3, 2))
par(mfrow=c(1, 1))

plot(1, type="n", main="",
     xlab="Event Time", xlim=c(0, 10), ylab="", yaxt='n', ylim=c(0, 1))

eps <- 0.05
ts <- c(-eps, t1, t2, t3, 10 + eps)
for (n in 1:4) {
  polygon(c(ts[n] + eps, ts[n + 1] - eps, ts[n + 1] - eps, ts[n] + eps), 
          c(0, 0, 1, 1), col = c_light, border = NA)
}

text(0, 0.9, cex=1.25, label="Left\nCensoring\nInterval",
     pos=4, col=c_dark)
text(3, 0.7, cex=1.25, label="Uncensored\nInterval",
     pos=4, col=c_dark)
text(6, 0.5, cex=1.25, label="Censored\nInterval",
     pos=4, col=c_dark)
text(8, 0.3, cex=1.25, label="Right\nCensoring\nInterval",
     pos=4, col=c_dark)

for (t in c(t1, t2, t3)) {
  abline(v=t, col="#DDDDDD", lty=2, lwd=3)
}


par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))

plot(1, type="n", main="",
     xlab="Event Time", xlim=c(0, 10), ylab="Interval Probability", ylim=c(0, 1))

eps <- 0.05
ts <- c(0, t1, t2, t3, 10)
for (n in 1:4) {
  p <- exp(-gamma * ts[n]) - exp(-gamma * ts[n + 1])
  polygon(c(ts[n] + eps, ts[n + 1] - eps, ts[n + 1] - eps, ts[n] + eps), 
          c(0, 0, p, p), col = c_light, border = NA)
}

for (t in c(t1, t2, t3)) {
  abline(v=t, col="#DDDDDD", lty=2, lwd=3)
}

# Simulate data using first method
simu_data <- list("N" = 350, "gamma" = gamma, 
                  "t1" = t1, "t2" = t2, "t3" = t3)

simu <- stan(file='stan_programs/simu_const_haz_cens1.stan', data=simu_data,
             iter=1, chains=1, seed=4838282, 
             algorithm="Fixed_param")

N1 <- extract(simu)$N1

N2 <- extract(simu)$N2
obs_times <- array(extract(simu)$obs_times[1,])

N3 <- extract(simu)$N3
N4 <- extract(simu)$N4

data <- list("t1" = t1, "t2" = t2, "t3" = t3,
             "N1" = N1, "N2" = N2, "N3" = N3, "N4" = N4,
             "obs_times" = obs_times)

# Visualize Data
buff_times <- c(rep(0.5 * t1, data$N1),
                data$obs_times,
                rep(0.5 * (t2 + t3), data$N3),
                rep(0.5 * (t4 + 10), data$N4))

par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))

plot_line_hist <- function(s, bins, xlab, main="", col="black") {
        B <- length(bins) - 1
        idx <- rep(1:B, each=2)
        x <- sapply(1:length(idx),
                    function(b) if(b %% 2 == 1) bins[idx[b]] else bins[idx[b] + 1])
        x <- c(0, x, 40)
        
        counts <- hist(s, breaks=bins, plot=FALSE)$counts
        y <- counts[idx]
        y <- c(0, y, 0)
        
        ymax <- max(y) + 1
        
        plot(x, y, type="l", main=main, col=col, lwd=2,
             xlab=xlab, xlim=c(min(bins), max(bins)),
             ylab="Counts", yaxt='n', ylim=c(0, ymax))
}

bins <- c(0, seq(2, 6, 0.5), 8, 10)
plot_line_hist(buff_times, bins, "Event Times")

for (t in c(t1, t2, t3)) {
        abline(v=t, col="#DDDDDD", lty=2, lwd=3)
}

par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))

M <- 100
xs <- seq(0, 10, 10 / (M - 1))
ys <- exp(-gamma * xs)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 1.01), ylab="Survival")

for (t in c(t1, t2, t3)) {
        abline(v=t, col="#DDDDDD", lty=2, lwd=3)
}

plot_eccdf(buff_times)


# Simulate data using second method
simu_data <- list("N" = 350, "gamma" = gamma, 
                  "t1" = t1, "t2" = t2, "t3" = t3)

simu <- stan(file='simu_const_haz_cens2.stan', data=simu_data,
             iter=1, chains=1, seed=4838282, 
             algorithm="Fixed_param")

N1 <- extract(simu)$N1

N2 <- extract(simu)$N2
obs_times <- array(extract(simu)$obs_times[1,])

N3 <- extract(simu)$N3
N4 <- extract(simu)$N4

data <- list("t1" = t1, "t2" = t2, "t3" = t3,
             "N1" = N1, "N2" = N2, "N3" = N3, "N4" = N4,
             "obs_times" = obs_times)

# Visualize Data
buff_times <- c(rep(0.5 * t1, data$N1),
                data$obs_times,
                rep(0.5 * (t2 + t3), data$N3),
                rep(0.5 * (t4 + 10), data$N4))

par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))

bins <- c(0, seq(2, 6, 0.5), 8, 10)
plot_line_hist(buff_times, bins, "Event Times")

for (t in c(t1, t2, t3)) {
        abline(v=t, col="#DDDDDD", lty=2, lwd=3)
}

par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))

M <- 100
xs <- seq(0, 10, 10 / (M - 1))
ys <- exp(-gamma * xs)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 1.01), ylab="Survival")

for (t in c(t1, t2, t3)) {
        abline(v=t, col="#DDDDDD", lty=2, lwd=3)
}

buff_times <- c(rep(0.5 * t1, N1),
                data$obs_times,
                rep(0.5 * (t2 + t3), N3),
                rep(0.5 * (t4 + 10), N4))

plot_eccdf(buff_times)


# Simulate data using first method conditonal on a fixed N2
simu_data <- list("N2" = 100, "gamma" = gamma, 
                  "t1" = t1, "t2" = t2, "t3" = t3)

simu <- stan(file='simu_const_haz_cens_cond1.stan', data=simu_data,
             iter=1, chains=1, seed=4838282, 
             algorithm="Fixed_param")

N1 <- extract(simu)$N1

obs_times <- array(extract(simu)$obs_times[1,])

N3 <- extract(simu)$N3
N4 <- extract(simu)$N4

data <- list("t1" = t1, "t2" = t2, "t3" = t3,
             "N1" = N1, "N2" = simu_data$N2, "N3" = N3, "N4" = N4,
             "obs_times" = obs_times)

# Visualize Data
buff_times <- c(rep(0.5 * t1, data$N1),
                data$obs_times,
                rep(0.5 * (t2 + t3), data$N3),
                rep(0.5 * (t4 + 10), data$N4))

par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))

bins <- c(0, seq(2, 6, 0.5), 8, 10)
plot_line_hist(buff_times, bins, "Event Times")

for (t in c(t1, t2, t3)) {
        abline(v=t, col="#DDDDDD", lty=2, lwd=3)
}

par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))

M <- 100
xs <- seq(0, 10, 10 / (M - 1))
ys <- exp(-gamma * xs)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 1.01), ylab="Survival")

for (t in c(t1, t2, t3)) {
        abline(v=t, col="#DDDDDD", lty=2, lwd=3)
}

buff_times <- c(rep(0.5 * t1, N1),
                data$obs_times,
                rep(0.5 * (t2 + t3), N3),
                rep(0.5 * (t4 + 10), N4))

plot_eccdf(buff_times)

# Simulate data using second method conditonal on a fixed N2
simu_data <- list("N2" = 100, "gamma" = gamma, 
                  "t1" = t1, "t2" = t2, "t3" = t3)

simu <- stan(file='simu_const_haz_cens_cond2.stan', data=simu_data,
             iter=1, chains=1, seed=4838282, 
             algorithm="Fixed_param")

N1 <- extract(simu)$N1

obs_times <- array(extract(simu)$obs_times[1,])

N3 <- extract(simu)$N3
N4 <- extract(simu)$N4

data <- list("t1" = t1, "t2" = t2, "t3" = t3,
             "N1" = N1, "N2" = simu_data$N2, "N3" = N3, "N4" = N4,
             "obs_times" = obs_times)

# Visualize Data
buff_times <- c(rep(0.5 * t1, data$N1),
                data$obs_times,
                rep(0.5 * (t2 + t3), data$N3),
                rep(0.5 * (t4 + 10), data$N4))

par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))

bins <- c(0, seq(2, 6, 0.5), 8, 10)
plot_line_hist(buff_times, bins, "Event Times")

for (t in c(t1, t2, t3)) {
  abline(v=t, col="#DDDDDD", lty=2, lwd=3)
}

par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))

M <- 100
xs <- seq(0, 10, 10 / (M - 1))
ys <- exp(-gamma * xs)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 1.01), ylab="Survival")

for (t in c(t1, t2, t3)) {
        abline(v=t, col="#DDDDDD", lty=2, lwd=3)
}

buff_times <- c(rep(0.5 * t1, N1),
                data$obs_times,
                rep(0.5 * (t2 + t3), N3),
                rep(0.5 * (t4 + 10), N4))

plot_eccdf(buff_times)

# Simulate final data for fitting
simu_data <- list("N" = 350, "gamma" = gamma, 
                  "t1" = t1, "t2" = t2, "t3" = t3)

simu <- stan(file='simu_const_haz_cens1.stan', data=simu_data,
             iter=1, chains=1, seed=4838282, 
             algorithm="Fixed_param")

N1 <- extract(simu)$N1[1]

N2 <- extract(simu)$N2[1]
obs_times <- array(extract(simu)$obs_times[1,])

N3 <- extract(simu)$N3[1]
N4 <- extract(simu)$N4[1]

data <- list("t1" = t1, "t2" = t2, "t3" = t3,
             "N1" = N1, "N2" = N2, "N3" = N3, "N4" = N4,
             "obs_times" = obs_times)

# Plot inferred survival behavior
data$N_display <- 100
data$display_times <- seq(0, 10, 10 / (data$N_display - 1))

data$t_min <- 0
data$t_max <- 10

# Fit with survival model
fit <- stan(file='stan_programs/fit_const_haz_cens.stan', data=data, 
            seed=4938483, refresh=0)

# False positives due to constant survival at t = 0
util$check_all_diagnostics(fit)

surv_samples = extract(fit)

par(mar = c(5, 2, 3, 2))
par(mfrow=c(1, 1))

plot_marginal("gamma", surv_samples, gamma, c(0, 1))

par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))

plot_marginal_quantiles(data$display_times, surv_samples$display_survival,
                        xlab="Time", xlim=c(0, 10),
                        ylab="Survival", ylim=c(0, 1))

buff_times <- c(rep(0.5 * t1, N1),
                data$obs_times,
                rep(0.5 * (t2 + t3), N3),
                rep(0.5 * (t4 + 10), N4))

for (t in c(t1, t2, t3)) {
        abline(v=t, col="#DDDDDD", lty=2, lwd=3)
}

plot_eccdf(buff_times)

plot_hist_retro <- function(data, pred_samples, name, breaks, xlim) {
  B <- length(breaks) - 1
  idx <- rep(1:B, each=2)
  xs <- sapply(1:length(idx),
               function(b) if(b %% 2 == 1) breaks[idx[b]] else breaks[idx[b] + 1])
        
  obs <- hist(data, breaks=breaks, plot=FALSE)$counts
  pad_obs <- do.call(cbind, lapply(idx, function(n) obs[n]))
        
  post_pred <- sapply(1:4000,
                      function(n) hist(pred_samples[n,], breaks=breaks, plot=FALSE)$counts)
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:B, function(b) quantile(post_pred[b,], probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9, n]))
        
  plot(1, type="n", main="Posterior Retrodictive Check",
       xlim=xlim, xlab=name,
       ylim=c(0, max(c(obs, cred[9,]))), ylab="")
        
  polygon(c(xs, rev(xs)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  for (b in 1:B)
    lines(xs[(2 * b - 1):(2 * b)], pad_cred[5,(2 * b - 1):(2 * b)], 
          col=c_dark, lwd=2)
        
  lines(xs, pad_obs, col="white", lty=1, lw=2.5)
  lines(xs, pad_obs, col="black", lty=1, lw=2)
}

plot_hist_retro_res <- function(data, pred_samples, name, breaks, xlim) {
  B <- length(breaks) - 1
  idx <- rep(1:B, each=2)
  xs <- sapply(1:length(idx),
               function(b) if(b %% 2 == 1) breaks[idx[b]] else breaks[idx[b] + 1])
  
  obs <- hist(data, breaks=breaks, plot=FALSE)$counts
    
  post_pred <- sapply(1:4000,
                      function(n) hist(pred_samples[n,], breaks=breaks, plot=FALSE)$counts 
                                - obs)
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:B, function(b) quantile(post_pred[b,], probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9, n]))
        
  plot(1, type="n", main="Residuals",
       xlim=xlim, xlab=name,
       ylim=c(min(cred[1,]), max(cred[9,])), ylab="Posterior Predictive - Observations")
        
  abline(h=0, col="#DDDDDD", lty=2, lwd=2)
        
  polygon(c(xs, rev(xs)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  for (b in 1:B)
    lines(xs[(2 * b - 1):(2 * b)], pad_cred[5,(2 * b - 1):(2 * b)], 
          col=c_dark, lwd=2)
}

par(mfrow=c(1, 2))

bins <- c(0, seq(2, 6, 0.5), 8, 10)

plot_hist_retro(buff_times, surv_samples$pred_obs_times, 
                "Event Times", bins, c(0, 10))
for (t in c(t1, t2, t3)) {
        abline(v=t, col="#DDDDDD", lty=2, lwd=2)
}

plot_hist_retro_res(buff_times, surv_samples$pred_obs_times, 
                    "Event Times", bins, c(0, 10))
for (t in c(t1, t2, t3)) {
        abline(v=t, col="#DDDDDD", lty=2, lwd=2)
}



par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))

plot_hist_retro(buff_times, surv_samples$pred_obs_times, 
                "Event Times", bins, c(0, 10))
for (t in c(t1, t2, t3)) {
        abline(v=t, col="#DDDDDD", lty=2, lwd=2)
}

text(-0.15, 150, cex=1.25, label="Left\nCensoring\nInterval",
     pos=4, col=c_dark)
text(2, 150, cex=1.25, label="Uncensored\nInterval",
     pos=4, col=c_dark)
text(5.95, 150, cex=1.25, label="Censored\nInterval",
     pos=4, col=c_dark)
text(8, 150, cex=1.25, label="Right\nCensoring\nInterval",
     pos=4, col=c_dark)


############################################################
# Prerequisite Events
############################################################

# Visualize True Behavior
gamma1 <- 0.2
gamma2 <- 0.5

par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 2))

M <- 100
xs <- seq(0, 35, 35 / (M - 1))
ys <- gamma1 * gamma2 * (exp(-gamma1 * xs) - exp(-gamma2 * xs)) / (gamma2 - gamma1)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 0.15), ylab="Event Density")

ys <- gamma2 / (gamma2 - gamma1) * exp(-gamma1 * xs) - gamma1 /(gamma2 - gamma1) * exp(-gamma2 * xs)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 1), ylab="Survival Probability")

ys <- exp(-gamma1 * xs)
lines(xs, ys, lwd=2, col=c_mid)

ys <- exp(-gamma2 * xs)
lines(xs, ys, lwd=2, col=c_light)

# Simulate Data
simu_data <- list("N" = 300, "gamma1" = gamma1, "gamma2" = gamma2)

simu <- stan(file='simu_const_haz_prereq.stan', data=simu_data,
             iter=1, chains=1, seed=4838282, 
             algorithm="Fixed_param")

obs_times <- array(extract(simu)$obs_times[1,])

data <- list("N" = simu_data$N,
             "obs_times" = obs_times)

# Visualize data
par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))

M <- 100
xs <- seq(0, 35, 35 / (M - 1))
ys <- exp(-gamma1 * xs)
plot(xs, ys, type="l", lwd=2, col=c_dark,
     ylim=c(0, 1), ylab="Survival")

ys <- exp(-gamma2 * xs)
lines(xs, ys, lwd=2, col=c_mid)

ys <- gamma2 / (gamma2 - gamma1) * exp(-gamma1 * xs) - gamma1 /(gamma2 - gamma1) * exp(-gamma2 * xs)
lines(xs, ys, lwd=2, col=c_light)

plot_eccdf(data$obs_times)

# Plot inferred survival behavior
data$N_display <- 100
data$display_times <- seq(0, 35, 35 / (data$N_display - 1))

# Fit with survival model
fit <- stan(file='fit_const_haz.stan', data=data, 
            seed=4938483, refresh=0)

# False positives due to constant survival at t = 0
util$check_all_diagnostics(fit)

surv_samples = extract(fit)

par(mar = c(5, 2, 3, 2))
par(mfrow=c(1, 1))

plot_marginal("gamma", surv_samples, NA, c(0, 0.5))

par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))
plot_marginal_quantiles(data$display_times, surv_samples$display_survival,
                        xlab="Time", xlim=c(0, 15),
                        ylab="Survival", ylim=c(0, 1))
plot_eccdf(data$obs_times)

plot_hist_retro <- function(data, pred_samples, name, min, max, delta,
                            xlim_override=NA) {
        
        breaks <- seq(min, max, delta)
        B <- length(breaks) - 1
        
        idx <- rep(1:B, each=2)
        xs <- sapply(1:length(idx),
                     function(b) if(b %% 2 == 1) breaks[idx[b]] else breaks[idx[b] + 1])
        
        obs <- hist(data, breaks=breaks, plot=FALSE)$counts
        pad_obs <- do.call(cbind, lapply(idx, function(n) obs[n]))
        
        post_pred <- sapply(1:4000,
                            function(n) hist(pred_samples[n,], breaks=breaks, plot=FALSE)$counts)
        probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
        cred <- sapply(1:B, function(b) quantile(post_pred[b,], probs=probs))
        pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9, n]))
        
        xlim <- c(min, max)
        if (!is.na(xlim_override))
                xlim <- xlim_override
        
        plot(1, type="n", main="Posterior Retrodictive Check",
             xlim=xlim, xlab=name,
             ylim=c(0, max(c(obs, cred[9,]))), ylab="")
        
        polygon(c(xs, rev(xs)), c(pad_cred[1,], rev(pad_cred[9,])),
                col = c_light, border = NA)
        polygon(c(xs, rev(xs)), c(pad_cred[2,], rev(pad_cred[8,])),
                col = c_light_highlight, border = NA)
        polygon(c(xs, rev(xs)), c(pad_cred[3,], rev(pad_cred[7,])),
                col = c_mid, border = NA)
        polygon(c(xs, rev(xs)), c(pad_cred[4,], rev(pad_cred[6,])),
                col = c_mid_highlight, border = NA)
        for (b in 1:B)
                lines(xs[(2 * b - 1):(2 * b)], pad_cred[5,(2 * b - 1):(2 * b)], 
                      col=c_dark, lwd=2)
        
        lines(xs, pad_obs, col="white", lty=1, lw=2.5)
        lines(xs, pad_obs, col="black", lty=1, lw=2)
        
        post_pred <- sapply(1:4000,
                            function(n) hist(pred_samples[n,], breaks=breaks, plot=FALSE)$counts 
                            - obs)
        probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
        cred <- sapply(1:B, function(b) quantile(post_pred[b,], probs=probs))
        pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9, n]))
        
        xlim <- c(min, max)
        if (!is.na(xlim_override))
                xlim <- xlim_override
        
        plot(1, type="n", main="Residuals",
             xlim=xlim, xlab=name,
             ylim=c(min(cred[1,]), max(cred[9,])), ylab="Posterior Predictive - Observations")
        
        abline(h=0, col="#DDDDDD", lty=2, lwd=2)
        
        polygon(c(xs, rev(xs)), c(pad_cred[1,], rev(pad_cred[9,])),
                col = c_light, border = NA)
        polygon(c(xs, rev(xs)), c(pad_cred[2,], rev(pad_cred[8,])),
                col = c_light_highlight, border = NA)
        polygon(c(xs, rev(xs)), c(pad_cred[3,], rev(pad_cred[7,])),
                col = c_mid, border = NA)
        polygon(c(xs, rev(xs)), c(pad_cred[4,], rev(pad_cred[6,])),
                col = c_mid_highlight, border = NA)
        for (b in 1:B)
                lines(xs[(2 * b - 1):(2 * b)], pad_cred[5,(2 * b - 1):(2 * b)], 
                      col=c_dark, lwd=2)
}

par(mfrow=c(1, 2))
plot_hist_retro(data$obs_times, surv_samples$pred_obs_times, 
                "Event Times", 0, 100, 2.5, c(0, 35))

# Fit with survival model with prerequisite events modeled with latent parameters
fit <- stan(file='fit_const_haz_prereq1.stan', data=data, 
            seed=4938483, refresh=0)

util$check_all_diagnostics(fit)

print(fit, pars=c("gamma1", "gamma2"))

surv_samples = extract(fit)

par(mar = c(5, 2, 3, 2))
par(mfrow=c(1, 2))

plot_marginal("gamma1", surv_samples, gamma1, c(0, 1))
plot_marginal("gamma2", surv_samples, gamma2, c(0, 1))

par(mfrow=c(1, 2))
plot_hist_retro(data$obs_times, surv_samples$pred_obs_times, 
                "Event Times", 0, 80, 2.5, c(0, 35))

# Fit with survival model
fit <- stan(file='fit_const_haz_prereq2.stan', data=data, 
            seed=4938483, refresh=0)

# False positives due to constant survival at t = 0
util$check_all_diagnostics(fit)

surv_samples = extract(fit)

par(mar = c(5, 2, 3, 2))
par(mfrow=c(1, 2))

plot_marginal("gamma1", surv_samples, gamma1, c(0, 1))
plot_marginal("gamma2", surv_samples, gamma2, c(0, 2))

par(mar = c(5, 4, 3, 2))
par(mfrow=c(1, 1))
plot_marginal_quantiles(data$display_times, surv_samples$display_survival,
                        xlab="Time", xlim=c(0, 30),
                        ylab="Survival", ylim=c(0, 1))
plot_eccdf(data$obs_times)

par(mfrow=c(1, 2))
plot_hist_retro(data$obs_times, surv_samples$pred_obs_times, 
                "Event Times", 0, 70, 2.5, c(0, 35))
