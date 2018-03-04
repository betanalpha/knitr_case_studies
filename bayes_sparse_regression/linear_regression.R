setwd('/Users/Betancourt/Documents/Research/Code/betanalpha/knitr_case_studies/bayes_sparse_regression')

############################################################
# Initial setup
############################################################

c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source('stan_utility.R')
source('plot_utility.R')

breaks <- function(mag, N) {
  return(mag*(2*(0:N)/N - 1))
}

############################################################
# Create data
############################################################

fit <- stan(file='generate_data.stan', iter=1,
            chains=1, seed=194838, algorithm="Fixed_param")

X <- extract(fit)$X[1,,]
y <- extract(fit)$y[1,]
N <- dim(X)[2]
M <- dim(X)[1]
beta_true <- extract(fit)$beta[1,]

stan_rdump(c("N", "M", "X", "y", "beta_true"), file="linear_regression.data.R")

beta_true <- extract(fit)$beta[1,]

############################################################
# Plot true slopes
############################################################

input_data <- read_rdump("linear_regression.data.R")

hist(input_data$beta_true, main="", col=c_dark, border=c_dark_highlight,
     xlab="True Slopes", yaxt='n', ylim=c(0, 75), ylab="", breaks=breaks(11, 100))

############################################################
# Linear regression with no priors on the slopes
############################################################

fit <- stan(file='linear_regression_unif.stan', data=input_data, seed=4938483)

# Check diagnostics
# Huge autocorrelations due to the non-identifiabilty
capture.output(check_n_eff(fit))[1:5]
capture.output(check_rhat(fit))[1:5]
check_div(fit)
check_treedepth(fit)
check_energy(fit)

# Extreme slopes due to non-identified likelihood
plot_post_quantiles(fit, input_data, "Uniform Prior")

# Reasonable recovery of truth
plot_residual_quantiles(fit, input_data, "Uniform Prior")

# Plot marginal posteriors for intercept and measurement variability
params = extract(fit)

par(mfrow=c(1, 2))

hist(params$alpha, main="", col=c_dark, border=c_dark_highlight,
xlab="alpha", yaxt='n', ylab="")
abline(v=3,col="white", lty=1, lwd=2.5)
abline(v=3,col="black", lty=1, lwd=2)

hist(params$sigma, main="", col=c_dark, border=c_dark_highlight,
xlab="sigma", yaxt='n', ylab="")
abline(v=1,col="white", lty=1, lwd=2.5)
abline(v=1,col="black", lty=1, lwd=2)

############################################################
# Linear regression with strongly regularizing priors on the slopes
############################################################

fit <- stan(file='linear_regression_narrow.stan', data=input_data, seed=4938483)

# Check diagnostics
# Much better fit with the stronger priors regularizing the non-identifiabilty
check_all_diagnostics(fit)

# Much more reasonable slopes
plot_post_quantiles(fit, input_data, "Narrow Prior")

# But the large slopes are not being recovered particularly well
# Our prior is regularizing the stronger slopes too much!
plot_residual_quantiles(fit, input_data, "Narrow Prior")

# Plot marginal posteriors for intercept and measurement variability
params = extract(fit)

par(mfrow=c(1, 2))

hist(params$alpha, main="", col=c_dark, border=c_dark_highlight,
xlab="alpha", yaxt='n', ylab="")
abline(v=3,col="white", lty=1, lwd=2.5)
abline(v=3,col="black", lty=1, lwd=2)

hist(params$sigma, main="", col=c_dark, border=c_dark_highlight,
xlab="sigma", yaxt='n', ylab="")
abline(v=1,col="white", lty=1, lwd=2.5)
abline(v=1,col="black", lty=1, lwd=2)

############################################################
# Linear regression with weakly regularizing priors on the slopes
############################################################

fit <- stan(file='linear_regression_wide.stan', data=input_data, seed=4938483)

# Check diagnostics
# Signs of nontrival autocorrelations due to the weak identifiability
check_all_diagnostics(fit)

# Reasonable slopes
plot_post_quantiles(fit, input_data, "Wide Prior")

# But the marginal posteriors are barely any more precise than the priors
# In particular the posteriors for the small slopes are extremely wide
plot_residual_quantiles(fit, input_data, "Wide Prior")

# Plot marginal posteriors for intercept and measurement variability
params = extract(fit)

par(mfrow=c(1, 2))

hist(params$alpha, main="", col=c_dark, border=c_dark_highlight,
xlab="alpha", yaxt='n', ylab="")
abline(v=3,col="white", lty=1, lwd=2.5)
abline(v=3,col="black", lty=1, lwd=2)

hist(params$sigma, main="", col=c_dark, border=c_dark_highlight,
xlab="sigma", yaxt='n', ylab="")
abline(v=1,col="white", lty=1, lwd=2.5)
abline(v=1,col="black", lty=1, lwd=2)

############################################################
# Linear regression with an L1 Prior
############################################################

fit <- stan(file='linear_regression_laplace.stan', data=input_data, seed=4938483,
            control=list(adapt_delta=0.99, max_treedepth=12))

# Check diagnostics
# Some indication of slow mixing over the measurement variability
check_n_eff(fit)
check_rhat(fit)
check_div(fit)
check_treedepth(fit, 12)
check_energy(fit)

# Much more reasonable slopes
plot_post_quantiles(fit, input_data, "Laplace Prior")

# Large slopes regularized too much, somewhat large uncertainties
plot_residual_quantiles(fit, input_data, "Laplace Prior")

# Plot marginal posteriors for intercept and measurement variability
params = extract(fit)

par(mfrow=c(1, 2))

hist(params$alpha, main="", col=c_dark, border=c_dark_highlight,
xlab="alpha", yaxt='n', ylab="")
abline(v=3,col="white", lty=1, lwd=2.5)
abline(v=3,col="black", lty=1, lwd=2)

hist(params$sigma, main="", col=c_dark, border=c_dark_highlight,
xlab="sigma", yaxt='n', ylab="")
abline(v=1,col="white", lty=1, lwd=2.5)
abline(v=1,col="black", lty=1, lwd=2)

############################################################
# Linear regression with a horseshoe prior
############################################################

fit <- stan(file='linear_regression_horseshoe.stan', data=input_data, seed=4938483,
            control=list(adapt_delta=0.99, max_treedepth=15))
# Only 74 divergences!
# Rhat(tau_tilde) is 1.17958342386435
# Rhat(sigma) is 1.1766476958234

fit <- stan(file='linear_regression_horseshoe2.stan', data=input_data, seed=4938483,
            control=list(adapt_delta=0.99, max_treedepth=15))
# 175 divergences
# Rhat(tau_a) = 1.34052825114609
# Rhat(sigma) = 1.24425464550976

# Check diagnostics
check_n_eff(fit)
check_rhat(fit)
check_div(fit)
check_treedepth(fit, 15)
check_energy(fit)

# Very narrow posteriors over the slopes
plot_post_quantiles(fit, input_data, "Horseshoe Prior")

# Not as consisent as we might like
plot_residual_quantiles(fit, input_data, "Horseshoe Prior")

# Plot marginal posteriors for intercept and measurement variability
params = extract(fit)

par(mfrow=c(1, 2))

hist(params$alpha, main="", col=c_dark, border=c_dark_highlight,
     xlab="alpha", yaxt='n', ylab="")
abline(v=3,col="white", lty=1, lwd=2.5)
abline(v=3,col="black", lty=1, lwd=2)

hist(params$sigma, main="", col=c_dark, border=c_dark_highlight,
     xlab="sigma", yaxt='n', ylab="")
abline(v=1,col="white", lty=1, lwd=2.5)
abline(v=1,col="black", lty=1, lwd=2)

############################################################
# Linear regression with the Finnish horseshoe prior
############################################################

fit <- stan(file='linear_regression_finnish_horseshoe.stan', data=input_data, seed=4938483,
            control=list(adapt_delta=0.99, max_treedepth=15))
# All diagnostics clear

fit <- stan(file='linear_regression_finnish_horseshoe2.stan', data=input_data, seed=4938483,
            control=list(adapt_delta=0.9, max_treedepth=15))


# Check diagnostics
check_n_eff(fit)
check_rhat(fit)
check_div(fit)
check_treedepth(fit, 15)
check_energy(fit)

# Extremely precise posteriors for both the small and large slopes
plot_post_quantiles(fit, input_data, "Finnish Horseshoe Prior")

# Both small and large slopes recovered accurately
plot_residual_quantiles(fit, input_data, "Finnish Horseshoe Prior")

# Plot marginal posteriors for intercept and measurement variability
params = extract(fit)

par(mfrow=c(1, 2))

hist(params$alpha, main="", col=c_dark, border=c_dark_highlight,
xlab="alpha", yaxt='n', ylab="")
abline(v=3,col="white", lty=1, lwd=2.5)
abline(v=3,col="black", lty=1, lwd=2)

hist(params$sigma, main="", col=c_dark, border=c_dark_highlight,
xlab="sigma", yaxt='n', ylab="")
abline(v=1,col="white", lty=1, lwd=2.5)
abline(v=1,col="black", lty=1, lwd=2)
