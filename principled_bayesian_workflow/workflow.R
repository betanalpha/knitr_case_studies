setwd('/Users/Betancourt/Documents/Research/Code/betanalpha/knitr_case_studies/principled_bayesian_workflow')

############################################################
#
# Initial Setup
#
############################################################

library(rstan)
rstan_options(auto_write = TRUE)

library(foreach)
library(doParallel)

util <- new.env()
source('stan_utility.R', local=util)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_mid_trans <- c("#B97C7C88")
c_mid_highlight_trans <- c("#A2505088")
c_dark_trans <- c("#8F272788")
c_dark_highlight_trans <- c("#7C000088")

############################################################
#
# Simulate Truth
#
############################################################

N <- 1000
simu_data <- list("N" = N)

fit <- stan(file='generative_truth.stan', data=simu_data,
            iter=1, warmup=0, chains=1, refresh=1,
            seed=4838282, algorithm="Fixed_param")

y <- extract(fit)$y[1,]

stan_rdump(c("N", "y"), file="workflow.data.R")

############################################################
#
# Workflow
#
############################################################

############################################################
# PRIOR TO OBSERVATION
############################################################

############################################################
# 1. Conceptual Analysis
############################################################

# We are working with a suite of detectors that each record
# the same source over a given time interval. The source
# strength and detector response are not expected to vary
# significantly in time.  Each detector is identical and
# returns discrete counts.

############################################################
# 2. Define Observations
############################################################

# Mathematically our observation takes the form of integer
# counts, y, for each of the N detectors.  In the Stan
# modeling lanugage this is specified as

writeLines(readLines("fit_data.stan", n=4))

############################################################
# 3. Identify Relevant Summary Statistics
############################################################

# There are N components in each observation, one for each
# detector.  We could analyze each component independently,
# but because we assume that the detectors are all identical
# we can analyze their comprehensive responses with a histogram
# of their counts.  In other words we consider the histogram
# of detector counts _as the summary statistic_!

# In this conceptual example assume that our conceptual
# domain expertise informs us that 25 counts in a detector
# would be an extreme but not impossible observation.

############################################################
# 4. Build a Generative Model
############################################################

# The constant source strength and detector responds suggests
# a Poisson observation model for each of the detectors with
# a single source strength, lambda.
#
# Our domain expertise that 25 counts is extreme suggests that
# we want our prior for lambda to keep most of its probability
# mass below lambda = 15, which corresponds to fluctations in
# the observations around 15 + 3 * sqrt(15) ~ 25.
#
# We achieve this with a half-normal prior with standard
# deviation = 6.44787 such that only 1% of the prior probability
# mass is above lambda = 15.

lambda <- seq(0, 20, 0.001)

plot(lambda, dnorm(lambda, 0, 6.44787), type="l", col=c_dark_highlight, lwd=2,
     xlab="lambda", ylab="Prior Density", yaxt='n')

lambda99 <- seq(0, 15, 0.001)
dens <- dnorm(lambda99, 0, 6.44787)
lambda99 <- c(lambda99, 15, 0)
dens <- c(dens, 0, 0)

polygon(lambda99, dens, col=c_dark, border=NA)

# This generative model is implemented in the Stan programs

writeLines(readLines("generative_ensemble.stan"))
writeLines(readLines("fit_data.stan"))

############################################################
# 5. Analyze the Generative Ensemble
############################################################

R <- 1000 # 1000 draws from the Bayesian joint distribution
N <- 100

############################################################
# 5a. Analyze the Prior Predictive Distribution
############################################################

simu_data <- list("N" = N)

fit <- stan(file='generative_ensemble.stan', data=simu_data,
            iter=R, warmup=0, chains=1, refresh=R,
            seed=4838282, algorithm="Fixed_param")

simu_lambdas <- extract(fit)$lambda
simu_ys <- extract(fit)$y

# Plot aggregated summary histogram for simulated observations
B <- 40
counts <- sapply(1:R, function(r) hist(simu_ys[r,], breaks=(0:(B + 1))-0.5, plot=FALSE)$counts)
probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:(B + 1), function(b) quantile(counts[b,], probs=probs))

idx <- rep(0:B, each=2)
x <- sapply(1:length(idx), function(b) if(b %% 2 == 0) idx[b] + 0.5 else idx[b] - 0.5)
pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n + 1]))

plot(1, type="n", main="Prior Predictive Distribution",
     xlim=c(-0.5, B + 0.5), xlab="y", ylim=c(0, max(cred[9,])), ylab="")

polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
        col = c_light, border = NA)
polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
        col = c_light_highlight, border = NA)
polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
        col = c_mid, border = NA)
polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
        col = c_mid_highlight, border = NA)
lines(x, pad_cred[5,], col=c_dark, lwd=2)

abline(v=25, col="white", lty=1, lw=2.5)
abline(v=25, col="black", lty=1, lw=2)

# We see a very small prior predictive probability above the
# extreme observation scale from our domain expertise

length(simu_ys[simu_ys > 25]) / length(simu_ys)

############################################################
# 5b. Fit the simulated observations and evaluate
############################################################

registerDoParallel(makeCluster(detectCores()))

simu_list <- t(data.matrix(data.frame(simu_lambdas, simu_ys)))

# Compile the posterior fit model
fit_model = stan_model(file='fit_data.stan')

ensemble_output <- foreach(simu=simu_list,
                           .combine='cbind') %dopar% {
  simu_lambda <- simu[1]
  simu_y <- simu[2:(N + 1)];

  # Fit the simulated observation
  input_data <- list("N" = N, "y" = simu_y)

  sink(file="/dev/null")
  library(rstan)
  fit <- sampling(fit_model, data=input_data, seed=4938483)
  sink()

  # Compute diagnostics
  util <- new.env()
  source('stan_utility.R', local=util)

  warning_code <- util$check_all_diagnostics(fit, quiet=TRUE)

  # Compute rank of prior draw with respect to thinned posterior draws
  sbc_rank <- sum(simu_lambda < extract(fit)$lambda[seq(1, 4000 - 8, 8)])

  # Compute posterior sensitivities
  s <- summary(fit, probs = c(), pars='lambda')$summary
  mean_lambda <- s[,1]
  sd_lambda <- s[,3]

  z_score <- abs((mean_lambda - simu_lambda) / sd_lambda)
  shrinkage <- 1 - sd_lambda * sd_lambda / 25.0

  c(warning_code, sbc_rank, z_score, shrinkage)
}

stopImplicitCluster()

# Check for fit diagnostics
warning_code <- ensemble_output[1,]
if (sum(warning_code) != 0) {
  print ("Some simulated posterior fits in the generative ensemble encountered problems!")
  for (r in 1:R) {
    if (warning_code[r] != 0) {
      print(sprintf('Replication %s of %s', r, R))
      util$parse_warning_code(warning_code[r])
      print(sprintf('Simulated lambda = %s', simu_lambdas[r]))
      print(" ")
    }
  }
} else {
  print ("No posterior fits in the generative ensemble encountered problems!")
}

# Check SBC histogram
sbc_rank <- ensemble_output[2,]
sbc_hist <- hist(sbc_rank, seq(0, 500, 25) - 0.5, col=c_dark, border=c_dark_highlight)
plot(sbc_hist, main="", xlab="Prior Rank", yaxt='n', ylab="")

low <- qbinom(0.005, R, 1 / 20)
mid <- qbinom(0.5, R, 1 / 20)
high <- qbinom(0.995, R, 1 / 20)
bar_x <- c(-10, 510, 500, 510, -10, 0, -10)
bar_y <- c(high, high, mid, low, low, mid, high)

polygon(bar_x, bar_y, col=c("#DDDDDD"), border=NA)
segments(x0=0, x1=500, y0=mid, y1=mid, col=c("#999999"), lwd=2)

plot(sbc_hist, col=c_dark, border=c_dark_highlight, add=T)

# Plot posterior sensitivities
z_score <- ensemble_output[3,]
shrinkage <- ensemble_output[4,]

plot(shrinkage, z_score, col=c("#8F272720"), lwd=2, pch=16, cex=0.8,
     xlim=c(0, 1), xlab="Posterior Shrinkage", ylim=c(0, 5), ylab="Posterior z-Score")

############################################################
# POSTERIOR TO OBSERVATION
############################################################

############################################################
# 6. Fit the observations and evaluate
############################################################

input_data <- read_rdump('workflow.data.R')
fit <- stan(file='fit_data_ppc.stan', data=input_data, seed=4938483)

# Check diagnostics
util$check_all_diagnostics(fit)

# Plot marginal posterior
params = extract(fit)

hist(params$lambda, main="", xlab="lambda", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

############################################################
# 7. Analyze the Posterior Predictive Distribution
############################################################

B <- 30

obs_counts <- hist(input_data$y, breaks=(0:(B + 1))-0.5, plot=FALSE)$counts

idx <- rep(0:B, each=2)
x <- sapply(1:length(idx), function(b) if(b %% 2 == 0) idx[b] + 0.5 else idx[b] - 0.5)
pad_obs <- do.call(cbind, lapply(idx, function(n) obs_counts[n + 1]))

counts <- sapply(1:4000, function(n) hist(params$y_ppc[n,], breaks=(0:(B + 1))-0.5, plot=FALSE)$counts)
probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:(B + 1), function(b) quantile(counts[b,], probs=probs))
pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n + 1]))

plot(1, type="n", main="Posterior Predictive Distribution",
     xlim=c(-0.5, B + 0.5), xlab="y",
     ylim=c(0, max(c(obs_counts, cred[9,]))), ylab="")

polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
        col = c_light, border = NA)
polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
        col = c_light_highlight, border = NA)
polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
        col = c_mid, border = NA)
polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
        col = c_mid_highlight, border = NA)
lines(x, pad_cred[5,], col=c_dark, lwd=2)

lines(x, pad_obs, col="white", lty=1, lw=2.5)
lines(x, pad_obs, col="black", lty=1, lw=2)

# The posterior predictive check indicates a serious
# excess of zeros above what we'd expect from a Poisson
# model.  Hence we want to expand our model to include
# zero-inflation and repeat the workflow.

############################################################
#
# Workflow Take 2
#
############################################################

############################################################
# PRIOR TO OBSERVATION
############################################################

############################################################
# 1. Conceptual Analysis
############################################################

# Our initial fit indicates that our detectors randomly
# malfunction and return zero counts.

############################################################
# 2. Define Observations
############################################################

# The observation space does not change.

writeLines(readLines("fit_data2.stan", n=4))

############################################################
# 3. Identify Relevant Summary Statistics
############################################################

# The same histogram summary from before will be productive
# with zero inflation.

############################################################
# 4. Build a Generative Model
############################################################

fit <- stan(file='prior_tune.stan', data=simu_data,
            iter=1, warmup=0, chains=1, refresh=R,
            seed=4838282, algorithm="Fixed_param")

# ZERO COMPONENT AND POISSON COMPONENT NON-IDENTIFIED IF
# THE COUNTS ARE VERY SMALL.  TO IDENTIFY THE MODEL WE
# NEED TO KEEP THE POISSON COMPONENT TO LAMBDA > 1.

# DRAW UPDATED LAMBDA AND THETA PRIORS!

# Our generative model expands to include zero-inflation
# by mixing the Poisson distribution with a Dirac
# distribution at zero.
#
# We have to be careful with this expansion, however,
# as it introduces a potential non-identifiability in
# our posterior; if the Poisson source strength is
# small then the Poisson distribution will be
# indistinguishible from the Dirac distribution at
# zero.  Consequently we want our prior for lambda
# to exclude this potential non-identifiabilty.  We
# do this with an inverse Gamma distribution that
# places only 1% probability below lambda = 1 and 1%
# above lambda = 15.

lambda <- seq(0, 20, 0.001)
dens <- lapply(lambda, function(l) dgamma(1 / l, 3.48681, rate=9.21604) * (1 / l**2))
plot(lambda, dens, type="l", col=c_dark_highlight, lwd=2,
     xlab="lambda", ylab="Prior Density", yaxt='n')

lambda98 <- seq(1, 15, 0.001)
dens <- lapply(lambda98, function(l) dgamma(1 / l, 3.48681, rate=9.21604) * (1 / l**2))
lambda98 <- c(lambda98, 15, 1)
dens <- c(dens, 0, 0)

polygon(lambda98, dens, col=c_dark, border=NA)

# We set a prior distribuiton for the zero mixture
# probability, theta, that puts only 1% probabilty
# mass below 0.1 and above 0.9 to ensure reasonable
# zero inflation.

theta <- seq(0, 1, 0.001)
dens <- dbeta(theta, 2.8663, 2.8663)
plot(theta, dens, type="l", col=c_dark_highlight, lwd=2,
     xlab="theta", ylab="Prior Density", yaxt='n')

theta98 <- seq(0.1, 0.9, 0.001)
dens <- dbeta(theta98, 2.8663, 2.8663)
theta98 <- c(theta98, 0.9, 0.1)
dens <- c(dens, 0, 0)

polygon(theta98, dens, col=c_dark, border=NA)

# This generative model is implemented in the Stan programs

writeLines(readLines("generative_ensemble2.stan"))
writeLines(readLines("fit_data2.stan"))

############################################################
# 5. Analyze the Generative Ensemble
############################################################

R <- 1000 # 1000 draws from the Bayesian joint distribution
N <- 100

############################################################
# 5a. Analyze the Prior Predictive Distribution
############################################################

simu_data <- list("N" = N)

fit <- stan(file='generative_ensemble2.stan', data=simu_data,
            iter=R, warmup=0, chains=1, refresh=R,
            seed=4838282, algorithm="Fixed_param")

simu_lambdas <- extract(fit)$lambda
simu_thetas <- extract(fit)$theta
simu_ys <- extract(fit)$y

# Plot aggregated summary histogram for simulated observations
hist(simu_ys, main="Prior Predictive Distribution", breaks=(0:70)-0.5,
     col=c_dark, border=c_dark_highlight,
     xlab="y", yaxt='n', ylab="")
abline(v=25, col=c_light, lty=1, lw=2)

# We see a very small prior predictive probability above the
# extreme observation scale from our domain expertise

length(simu_ys[simu_ys > 25]) / length(simu_ys)

############################################################
# 5b. Fit the simulated observations and evaluate
############################################################

registerDoParallel(makeCluster(detectCores()))

simu_list <- t(data.matrix(data.frame(simu_lambdas, simu_thetas, simu_ys)))

# Compile the posterior fit model
fit_model = stan_model(file='fit_data2.stan')

ensemble_output <- foreach(simu=simu_list,
                           .combine='cbind') %dopar% {
  simu_lambda <- simu[1]
  simu_theta <- simu[2]
  simu_y <- simu[3:(N + 2)];

  # Fit the simulated observation
  input_data <- list("N" = N, "y" = simu_y)

  sink(file="/dev/null")
  library(rstan)
  fit <- sampling(fit_model, data=input_data, seed=4938483)
  sink()

  # Compute diagnostics
  util <- new.env()
  source('stan_utility.R', local=util)

  warning_code <- util$check_all_diagnostics(fit, quiet=TRUE)

  # Compute rank of prior draw with respect to thinned posterior draws
  sbc_rank_lambda <- sum(simu_lambda < extract(fit)$lambda[seq(1, 4000 - 8, 8)])
  sbc_rank_theta <- sum(simu_theta < extract(fit)$theta[seq(1, 4000 - 8, 8)])

  # Compute posterior sensitivities
  s <- summary(fit, probs = c(), pars='lambda')$summary
  mean_lambda <- s[,1]
  sd_lambda <- s[,3]

  z_score_lambda <- abs((mean_lambda - simu_lambda) / sd_lambda)
  shrinkage_lambda <- 1 - sd_lambda * sd_lambda / 25.0

  s <- summary(fit, probs = c(), pars='theta')$summary
  mean_theta <- s[,1]
  sd_theta <- s[,3]

  z_score_theta <- abs((mean_theta - simu_theta) / sd_theta)
  shrinkage_theta <- 1 - sd_theta * sd_theta / 0.288

  c(warning_code,
    sbc_rank_lambda, z_score_lambda, shrinkage_lambda,
    sbc_rank_theta, z_score_theta, shrinkage_theta)
}

stopImplicitCluster()

# Check for fit diagnostics
warning_code <- ensemble_output[1,]
if (sum(warning_code) != 0) {
  print ("Some simulated posterior fits in the generative ensemble encountered problems!")
  for (r in 1:R) {
    if (warning_code[r] != 0) {
      print(sprintf('Replication %s of %s', r, R))
      util$parse_warning_code(warning_code[r])
      print(sprintf('Simulated lambda = %s', simu_lambdas[r]))
      print(sprintf('Simulated theta = %s', simu_thetas[r]))
      print(" ")
    }
  }
} else {
  print ("No posterior fits in the generative ensemble encountered problems!")
}

# Check SBC histogram
sbc_rank <- ensemble_output[2,]
sbc_hist <- hist(sbc_rank, seq(0, 500, 25) - 0.5, col=c_dark, border=c_dark_highlight)
plot(sbc_hist, main="", xlab="Prior Rank (Lambda)", yaxt='n', ylab="")

low <- qbinom(0.005, R, 1 / 20)
mid <- qbinom(0.5, R, 1 / 20)
high <- qbinom(0.995, R, 1 / 20)
bar_x <- c(-10, 510, 500, 510, -10, 0, -10)
bar_y <- c(high, high, mid, low, low, mid, high)

polygon(bar_x, bar_y, col=c("#DDDDDD"), border=NA)
segments(x0=0, x1=500, y0=mid, y1=mid, col=c("#999999"), lwd=2)

plot(sbc_hist, col=c_dark, border=c_dark_highlight, add=T)

sbc_rank <- ensemble_output[5,]
sbc_hist <- hist(sbc_rank, seq(0, 500, 25) - 0.5, col=c_dark, border=c_dark_highlight)
plot(sbc_hist, main="", xlab="Prior Rank (Lambda)", yaxt='n', ylab="")

mean <- length(sbc_rank) / (500 / 25)
low <- mean - 3 * sqrt(mean)
high <- mean + 3 * sqrt(mean)
bar_x <- c(-10, 510, 500, 510, -10, 0, -10)
bar_y <- c(high, high, mean, low, low, mean, high)

polygon(bar_x, bar_y, col=c("#DDDDDD"), border=NA)
segments(x0=0, x1=500, y0=mean, y1=mean, col=c("#999999"), lwd=2)

plot(sbc_hist, col=c_dark, border=c_dark_highlight, add=T)

# Plot posterior sensitivities
z_score <- ensemble_output[3,]
shrinkage <- ensemble_output[4,]

plot(shrinkage, z_score, col=c("#8F272720"), lwd=2, pch=16, cex=0.8, main="Lambda",
     xlim=c(0, 1), xlab="Posterior Shrinkage", ylim=c(0, 5), ylab="Posterior z-Score")

z_score <- ensemble_output[6,]
shrinkage <- ensemble_output[7,]

plot(shrinkage, z_score, col=c("#8F272720"), lwd=2, pch=16, cex=0.8, main="Theta",
     xlim=c(0, 1), xlab="Posterior Shrinkage", ylim=c(0, 5), ylab="Posterior z-Score")

############################################################
# POSTERIOR TO OBSERVATION
############################################################

############################################################
# 6. Fit the observations and evaluate
############################################################

input_data <- read_rdump('workflow.data.R')
fit <- stan(file='fit_data2_ppc.stan', data=input_data, seed=4938483)

# Check diagnostics
util$check_all_diagnostics(fit)

# Plot marginal posterior
params = extract(fit)

par(mfrow=c(2, 1))

hist(params$lambda, main="", xlab="lambda", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(params$theta, main="", xlab="theta", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

############################################################
# 7. Analyze the Posterior Predictive Distribution
############################################################

y_ppc = params$y_ppc
dim(y_ppc) <- NULL

p1 <- hist(y_ppc, breaks=(0:25))
p1$counts = p1$counts / sum(p1$counts)
p2 <- hist(input_data$y, breaks=(0:25))
p2$counts = p2$counts / sum(p2$counts)

plot(p2, main="Posterior Predictive Check", xlab="y", yaxt='n', ylab="",
col=c_dark_trans, border=c_dark_highlight_trans)
plot(p1, col=c_mid_trans, border=c_mid_highlight_trans, add=T)

legend("topright", c("Observed Data", "Posterior Predictive Distribution"),
fill=c(c_dark_trans, c_mid_trans), bty="n")

# Looks much better, but we can still see that the tail of
# posterior predictive distribution goes significantly past
# the observe data indicating that the detectors are in
# fact censored and do not report any values if the true
# counts are above y = 14.  We'll need to expand our model
# again and repeat the workflow.

############################################################
#
# Workflow Take 3
#
############################################################

############################################################
# PRIOR TO OBSERVATION
############################################################

############################################################
# 1. Conceptual Analysis
############################################################

# Our second fit indicates that our detectors not only
# report zeroes occassionaly but also fail to report
# any counts above y = 14 at all.  In other words we
# actually have a much larger suite of detectors and
# record only the first 100 that return a value below
# y = 14.

############################################################
# 2. Define Observations
############################################################

# The observation space does not change.

writeLines(readLines("fit_data3.stan", n=4))

############################################################
# 3. Identify Relevant Summary Statistics
############################################################

# The same histogram summary from before will be productive
# with zero inflation and censoring.

############################################################
# 4. Build a Generative Model
############################################################

# Our generative model expands to include a censored
# Poisson distribution in addition to zero inflation.

writeLines(readLines("generative_ensemble3.stan"))
writeLines(readLines("fit_data3.stan"))

############################################################
# 5. Analyze the Generative Ensemble
############################################################

R <- 1000 # 1000 draws from the Bayesian joint distribution
N <- 100

############################################################
# 5a. Analyze the Prior Predictive Distribution
############################################################

simu_data <- list("N" = N)

fit <- stan(file='generative_ensemble3.stan', data=simu_data,
            iter=R, warmup=0, chains=1, refresh=R,
            seed=4838282, algorithm="Fixed_param")

simu_lambdas <- extract(fit)$lambda
simu_thetas <- extract(fit)$theta
simu_ys <- extract(fit)$y

# Plot aggregated summary histogram for simulated observations
hist(simu_ys, main="Prior Predictive Distribution", breaks=(0:26)-0.5,
     col=c_dark, border=c_dark_highlight,
     xlab="y", yaxt='n', ylab="")
abline(v=25, col=c_light, lty=1, lw=2)

# With the censoring we now have exactly prior predictive
# probability above the extreme observation scale from our
# domain expertise

length(simu_ys[simu_ys > 25]) / length(simu_ys)

############################################################
# 5b. Fit the simulated observations and evaluate
############################################################

registerDoParallel(makeCluster(detectCores()))

simu_list <- t(data.matrix(data.frame(simu_lambdas, simu_thetas, simu_ys)))

# Compile the posterior fit model
fit_model = stan_model(file='fit_data3.stan')

ensemble_output <- foreach(simu=simu_list,
                           .combine='cbind') %dopar% {
  simu_lambda <- simu[1]
  simu_theta <- simu[2]
  simu_y <- simu[3:(N + 2)];

  # Fit the simulated observation
  input_data <- list("N" = N, "y" = simu_y)

  sink(file="/dev/null")
  library(rstan)
  fit <- sampling(fit_model, data=input_data, seed=4938483)
  sink()

  # Compute diagnostics
  util <- new.env()
  source('stan_utility.R', local=util)

  warning_code <- util$check_all_diagnostics(fit, quiet=TRUE)

  # Compute rank of prior draw with respect to thinned posterior draws
  sbc_rank_lambda <- sum(simu_lambda < extract(fit)$lambda[seq(1, 4000 - 8, 8)])
  sbc_rank_theta <- sum(simu_theta < extract(fit)$theta[seq(1, 4000 - 8, 8)])

  # Compute posterior sensitivities
  s <- summary(fit, probs = c(), pars='lambda')$summary
  mean_lambda <- s[,1]
  sd_lambda <- s[,3]

  z_score_lambda <- abs((mean_lambda - simu_lambda) / sd_lambda)
  shrinkage_lambda <- 1 - sd_lambda * sd_lambda / 25.0

  s <- summary(fit, probs = c(), pars='theta')$summary
  mean_theta <- s[,1]
  sd_theta <- s[,3]

  z_score_theta <- abs((mean_theta - simu_theta) / sd_theta)
  shrinkage_theta <- 1 - sd_theta * sd_theta / 0.288

  c(warning_code,
  sbc_rank_lambda, z_score_lambda, shrinkage_lambda,
  sbc_rank_theta, z_score_theta, shrinkage_theta)
}

stopImplicitCluster()

# Check for fit diagnostics
warning_code <- ensemble_output[1,]
if (sum(warning_code) != 0) {
  print ("Some simulated posterior fits in the generative ensemble encountered problems!")
  for (r in 1:R) {
    if (warning_code[r] != 0) {
      print(sprintf('Replication %s of %s', r, R))
      util$parse_warning_code(warning_code[r])
      print(sprintf('Simulated lambda = %s', simu_lambdas[r]))
      print(sprintf('Simulated theta = %s', simu_thetas[r]))
      print(" ")
    }
  }
} else {
  print ("No posterior fits in the generative ensemble encountered problems!")
}

# Check SBC histogram
sbc_rank <- ensemble_output[2,]
sbc_hist <- hist(sbc_rank, seq(0, 500, 25) - 0.5, col=c_dark, border=c_dark_highlight)
plot(sbc_hist, main="", xlab="Prior Rank (Lambda)", yaxt='n', ylab="")

low <- qbinom(0.005, R, 1 / 20)
mid <- qbinom(0.5, R, 1 / 20)
high <- qbinom(0.995, R, 1 / 20)
bar_x <- c(-10, 510, 500, 510, -10, 0, -10)
bar_y <- c(high, high, mid, low, low, mid, high)

polygon(bar_x, bar_y, col=c("#DDDDDD"), border=NA)
segments(x0=0, x1=500, y0=mid, y1=mid, col=c("#999999"), lwd=2)

plot(sbc_hist, col=c_dark, border=c_dark_highlight, add=T)

sbc_rank <- ensemble_output[5,]
sbc_hist <- hist(sbc_rank, seq(0, 500, 25) - 0.5, col=c_dark, border=c_dark_highlight)
plot(sbc_hist, main="", xlab="Prior Rank (Lambda)", yaxt='n', ylab="")

mean <- length(sbc_rank) / (500 / 25)
low <- mean - 3 * sqrt(mean)
high <- mean + 3 * sqrt(mean)
bar_x <- c(-10, 510, 500, 510, -10, 0, -10)
bar_y <- c(high, high, mean, low, low, mean, high)

polygon(bar_x, bar_y, col=c("#DDDDDD"), border=NA)
segments(x0=0, x1=500, y0=mean, y1=mean, col=c("#999999"), lwd=2)

plot(sbc_hist, col=c_dark, border=c_dark_highlight, add=T)

# Plot posterior sensitivities
z_score <- ensemble_output[3,]
shrinkage <- ensemble_output[4,]

plot(shrinkage, z_score, col=c("#8F272720"), lwd=2, pch=16, cex=0.8, main="Lambda",
     xlim=c(0, 1), xlab="Posterior Shrinkage", ylim=c(0, 5), ylab="Posterior z-Score")

z_score <- ensemble_output[6,]
shrinkage <- ensemble_output[7,]

plot(shrinkage, z_score, col=c("#8F272720"), lwd=2, pch=16, cex=0.8, main="Theta",
    xlim=c(0, 1), xlab="Posterior Shrinkage", ylim=c(0, 5), ylab="Posterior z-Score")

############################################################
# POSTERIOR TO OBSERVATION
############################################################

############################################################
# 6. Fit the observations and evaluate
############################################################

input_data <- read_rdump('workflow.data.R')
fit <- stan(file='fit_data3_ppc.stan', data=input_data, seed=4938483)

# Check diagnostics
util$check_all_diagnostics(fit)

# Plot marginal posterior
params = extract(fit)

par(mfrow=c(2, 1))

hist(params$lambda, main="", xlab="lambda", yaxt='n', ylab="",
col=c_dark, border=c_dark_highlight)

hist(params$theta, main="", xlab="theta", yaxt='n', ylab="",
col=c_dark, border=c_dark_highlight)

############################################################
# 7. Analyze the Posterior Predictive Distribution
############################################################

y_ppc = params$y_ppc
dim(y_ppc) <- NULL

p1 <- hist(y_ppc, breaks=(0:25))
p1$counts = p1$counts / sum(p1$counts)
p2 <- hist(input_data$y, breaks=(0:25))
p2$counts = p2$counts / sum(p2$counts)

plot(p2, main="Posterior Predictive Check", xlab="y", yaxt='n', ylab="",
col=c_dark_trans, border=c_dark_highlight_trans)
plot(p1, col=c_mid_trans, border=c_mid_highlight_trans, add=T)

legend("topright", c("Observed Data", "Posterior Predictive Distribution"),
fill=c(c_dark_trans, c_mid_trans), bty="n")

# We finally have a model that captures the intricacies
# or our observed data!
