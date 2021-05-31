############################################################
# Initial setup
############################################################

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
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 5))

############################################################
# Graphical Demonstrations
############################################################

delta <- 0.025
xs <- seq(-15, 15, delta)

density_comps <- data.frame(xs)

# Wide Likelihood Function
sigma <- 4
y <- 2

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

density_comps['wide_likelihoods'] = likelihoods

# Narrow Likelihood Function Concentrating At Small Values
sigma <- 1
y <- 1.5

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

density_comps['narrow_small_likelihoods'] = likelihoods

# Narrow Likelihood Function Concentrating At Large Values
sigma <- 1
y <- 10.25

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

density_comps['narrow_large_likelihoods'] = likelihoods

# Shrinkage
sigma <- 1

N <- 500
likelihood_peaks <- seq(0, 8, 8 / (N - 1))
unif_posterior_means <- rep(0, N)
unif_posterior_sds <- rep(0, N)

# Uniform Prior
for (n in 1:N) {
  unnorm_likelihoods <- exp(- 0.5 * ( (xs - likelihood_peaks[n]) / sigma)**2)
  
  unnorm_posterior_densities <- unnorm_likelihoods
  
  norm <- sum(unnorm_posterior_densities) * delta
  posterior_densities <- unnorm_posterior_densities / norm
  
  unif_posterior_means[n] <- sum(xs * posterior_densities) * delta
  unif_posterior_sds[n] <- sqrt(sum( (xs - unif_posterior_means[n])**2 * posterior_densities) * delta)
} 

shrinkage <- data.frame(likelihood_peaks, unif_posterior_means, unif_posterior_sds)

############################################################
# Normal
############################################################

par(mfrow=c(1, 1))

# Prior
delta <- 0.025
xs <- seq(-15, 15, delta)

beta <- 1

prior_densities <- exp(- 0.5 * (xs / beta)**2) / sqrt(2 * pi * beta**2)

density_comps['normal_prior'] = prior_densities

# Wide Likelihood Function
sigma <- 4
y <- 5

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

unnorm_posterior_densities <- unnorm_likelihoods * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, likelihoods, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

density_comps['normal_posterior_wide'] = posterior_densities

# Narrow Likelihood Function Concentrating At Small Values
sigma <- 1
y <- 1.5

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

unnorm_posterior_densities <- unnorm_likelihoods * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, likelihoods, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

density_comps['normal_posterior_narrow_small'] = posterior_densities

# Narrow Likelihood Function Concentrating At Large Values
sigma <- 1
y <- 10.25

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

unnorm_posterior_densities <- unnorm_likelihoods * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, likelihoods, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

density_comps['normal_posterior_narrow_large'] = posterior_densities

# Shrinkage
delta <- 0.005
xs <- seq(-12, 12, delta)
N_x <- length(xs)

beta <- 1

prior_densities <- exp(-0.5 * (xs / beta)**2) / sqrt(2 * pi * beta**2) 

sigma <- 1
N <- 500
posterior_means <- rep(0, N)
posterior_sds <- rep(0, N)

for (n in 1:N) {
  unnorm_likelihoods <- exp(- 0.5 * ( (xs - likelihood_peaks[n]) / sigma)**2)
  
  unnorm_posterior_densities <- unnorm_likelihoods * prior_densities
  
  norm <- sum(unnorm_posterior_densities) * delta
  posterior_densities <- unnorm_posterior_densities / norm
  
  posterior_means[n] <- sum(xs * posterior_densities) * delta
  posterior_sds[n] <- sqrt(sum( (xs - posterior_means[n])**2 * posterior_densities) * delta)
}

plot(likelihood_peaks, posterior_means / unif_posterior_means, 
     type="l", lwd=2, col=c_dark, xlim=c(5 * delta, 8), ylim=c(0, 1.1))
lines(c(0, 8), c(1, 1), lwd=2, col=c_light)
lines(c(beta, beta), c(0, 1.1), lwd=2, col=c_light)

plot(likelihood_peaks, posterior_sds / unif_posterior_sds, 
     type="l", lwd=2, col=c_dark, xlim=c(0, 8), ylim=c(0, 1.1))
lines(c(0, 8), c(1, 1), lwd=2, col=c_light)
lines(c(beta, beta), c(0, 1.1), lwd=2, col=c_light)

shrinkage['normal_means'] <- posterior_means
shrinkage['normal_sds'] <- posterior_sds

############################################################
# Normal Mixture
############################################################

par(mfrow=c(1, 1))

# Prior
delta <- 0.025
xs <- seq(-15, 15, delta)

beta1 <- 1
beta2 <- 10
gamma <- 0.75
  
prior_densities <- gamma * exp(-0.5 * (xs / beta1)**2) / sqrt(2 * pi * beta1**2) +
                   (1 - gamma) * exp(-0.5 * (xs / beta2)**2) / sqrt(2 * pi * beta2**2) 

density_comps['normal_mixture_prior'] = prior_densities

# Wide Likelihood Function
sigma <- 4
y <- 5

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

unnorm_posterior_densities <- unnorm_likelihoods * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, likelihoods, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

density_comps['normal_mixture_posterior_wide'] = posterior_densities

# Narrow Likelihood Function Concentrating At Small Values
sigma <- 1
y <- 1.5

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

unnorm_posterior_densities <- unnorm_likelihoods * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, likelihoods, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

density_comps['normal_mixture_posterior_narrow_small'] = posterior_densities

# Narrow Likelihood Function Concentrating At Large Values
sigma <- 1
y <- 10.25

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

unnorm_posterior_densities <- unnorm_likelihoods * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, likelihoods, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

density_comps['normal_mixture_posterior_narrow_large'] = posterior_densities

# Shrinkage
delta <- 0.005
xs <- seq(-15, 15, delta)
N_x <- length(xs)

beta <- 1

prior_densities <- gamma * exp(-0.5 * (xs / beta1)**2) / sqrt(2 * pi * beta1**2) +
                   (1 - gamma) * exp(-0.5 * (xs / beta2)**2) / sqrt(2 * pi * beta2**2) 

sigma <- 1

N <- 500
posterior_means <- rep(0, N)
posterior_sds <- rep(0, N)

for (n in 1:N) {
  unnorm_likelihoods <- exp(- 0.5 * ( (xs - likelihood_peaks[n]) / sigma)**2)
  
  unnorm_posterior_densities <- unnorm_likelihoods * prior_densities
  
  norm <- sum(unnorm_posterior_densities) * delta
  posterior_densities <- unnorm_posterior_densities / norm
  
  posterior_means[n] <- sum(xs * posterior_densities) * delta
  posterior_sds[n] <- sqrt(sum( (xs - posterior_means[n])**2 * posterior_densities) * delta)
}

plot(likelihood_peaks, posterior_means / unif_posterior_means, 
     type="l", lwd=2, col=c_dark, xlim=c(5 * delta, 8), ylim=c(0, 1.1))
lines(c(0, 8), c(1, 1), lwd=2, col=c_light)
lines(c(beta1, beta1), c(0, 1.1), lwd=2, col=c_light)
lines(c(beta2, beta2), c(0, 1.1), lwd=2, col=c_light)

plot(likelihood_peaks, posterior_sds / unif_posterior_sds, 
     type="l", lwd=2, col=c_dark, xlim=c(0, 8), ylim=c(0, 1.35))
lines(c(0, 8), c(1, 1), lwd=2, col=c_light)
lines(c(beta, beta), c(0, 1.1), lwd=2, col=c_light)
lines(c(beta2, beta2), c(0, 1.1), lwd=2, col=c_light)

shrinkage['normal_mixture_means'] <- posterior_means
shrinkage['normal_mixture_sds'] <- posterior_sds

############################################################
# Laplace
############################################################

par(mfrow=c(1, 1))

# Prior
delta <- 0.025
xs <- seq(-15, 15, delta)

beta <- 1

prior_densities <- exp( - abs(xs) / beta) / (2 * beta) 

density_comps['laplace_prior'] = prior_densities

# Wide Likelihood Function
sigma <- 4
y <- 5

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

unnorm_posterior_densities <- unnorm_likelihoods * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, likelihoods, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

density_comps['laplace_posterior_wide'] = posterior_densities

# Narrow Likelihood Function Concentrating At Small Values
sigma <- 1
y <- 1.5

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

unnorm_posterior_densities <- unnorm_likelihoods * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, likelihoods, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

density_comps['laplace_posterior_narrow_small'] = posterior_densities

# Narrow Likelihood Function Concentrating At Large Values
sigma <- 1
y <- 10.25

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

unnorm_posterior_densities <- unnorm_likelihoods * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, likelihoods, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

density_comps['laplace_posterior_narrow_large'] = posterior_densities

# Shrinkage
delta <- 0.005
xs <- seq(-15, 15, delta)
N_x <- length(xs)

beta <- 1

prior_densities <- exp( - abs(xs) / beta) / (2 * beta) 

sigma <- 1

N <- 500
posterior_means <- rep(0, N)
posterior_sds <- rep(0, N)

for (n in 1:N) {
  unnorm_likelihoods <- exp(- 0.5 * ( (xs - likelihood_peaks[n]) / sigma)**2)
  
  unnorm_posterior_densities <- unnorm_likelihoods * prior_densities
  
  norm <- sum(unnorm_posterior_densities) * delta
  posterior_densities <- unnorm_posterior_densities / norm
  
  posterior_means[n] <- sum(xs * posterior_densities) * delta
  posterior_sds[n] <- sqrt(sum( (xs - posterior_means[n])**2 * posterior_densities) * delta)
}

plot(likelihood_peaks, posterior_means / unif_posterior_means, 
     type="l", lwd=2, col=c_dark, xlim=c(5 * delta, 8), ylim=c(0, 1.1))
lines(c(0, 8), c(1, 1), lwd=2, col=c_light)
lines(c(beta, beta), c(0, 1.1), lwd=2, col=c_light)

plot(likelihood_peaks, posterior_sds / unif_posterior_sds, 
     type="l", lwd=2, col=c_dark, xlim=c(0, 8), ylim=c(0, 1.1))
lines(c(0, 8), c(1, 1), lwd=2, col=c_light)
lines(c(beta, beta), c(0, 1.1), lwd=2, col=c_light)

shrinkage['laplace_means'] <- posterior_means
shrinkage['laplace_sds'] <- posterior_sds

############################################################
# Cauchy
############################################################

# Prior
delta <- 0.025
xs <- seq(-15, 15, delta)

beta <- 1

prior_densities <- (beta**2 / (xs**2 + beta**2)) / (pi * beta) 

density_comps['cauchy_prior'] = prior_densities

# Wide Likelihood Function
sigma <- 4
y <- 5

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

unnorm_posterior_densities <- unnorm_likelihoods * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, likelihoods, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

density_comps['cauchy_posterior_wide'] = posterior_densities

# Narrow Likelihood Function Concentrating At Small Values
sigma <- 1
y <- 1.5

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

unnorm_posterior_densities <- unnorm_likelihoods * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, likelihoods, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

density_comps['cauchy_posterior_narrow_small'] = posterior_densities

# Narrow Likelihood Function Concentrating At Large Values
sigma <- 1
y <- 10.25

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

unnorm_posterior_densities <- unnorm_likelihoods * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, likelihoods, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

density_comps['cauchy_posterior_narrow_large'] = posterior_densities

# Shrinkage
delta <- 0.005
xs <- seq(-15, 15, delta)
N_x <- length(xs)

beta <- 1

prior_densities <- (beta**2 / (xs**2 + beta**2)) / (pi * beta) 

sigma <- 1

N <- 500
posterior_means <- rep(0, N)
posterior_sds <- rep(0, N)

for (n in 1:N) {
  unnorm_likelihoods <- exp(- 0.5 * ( (xs - likelihood_peaks[n]) / sigma)**2)
  
  unnorm_posterior_densities <- unnorm_likelihoods * prior_densities
  
  norm <- sum(unnorm_posterior_densities) * delta
  posterior_densities <- unnorm_posterior_densities / norm
  
  posterior_means[n] <- sum(xs * posterior_densities) * delta
  posterior_sds[n] <- sqrt(sum( (xs - posterior_means[n])**2 * posterior_densities) * delta)
}

plot(likelihood_peaks, posterior_means / unif_posterior_means, 
     type="l", lwd=2, col=c_dark, xlim=c(5 * delta, 8), ylim=c(0, 1.1))
lines(c(0, 8), c(1, 1), lwd=2, col=c_light)
lines(c(beta, beta), c(0, 1.1), lwd=2, col=c_light)

plot(likelihood_peaks, posterior_sds / unif_posterior_sds, 
     type="l", lwd=2, col=c_dark, xlim=c(0, 8), ylim=c(0, 1.1))
lines(c(0, 8), c(1, 1), lwd=2, col=c_light)
lines(c(beta, beta), c(0, 1.1), lwd=2, col=c_light)

shrinkage['cauchy_means'] <- posterior_means
shrinkage['cauchy_sds'] <- posterior_sds

############################################################
# Save Output
############################################################

write.csv(density_comps, "density_comps.csv", row.names=FALSE, quote=FALSE)
write.csv(shrinkage, "shrinkage.csv", row.names=FALSE, quote=FALSE)


############################################################
#
# Outer Extent Cumulant Regularization
#
############################################################

delta <- 0.01
xs <- seq(-300, 300, delta)
N_x <- length(xs)

sigma <- 20

N <- 500
likelihood_peaks <- seq(0, 200, 200 / (N - 1))
unif_posterior_means <- rep(0, N)
unif_posterior_sds <- rep(0, N)

# Uniform Prior
for (n in 1:N) {
  unnorm_likelihoods <- exp(- 0.5 * ( (xs - likelihood_peaks[n]) / sigma)**2)
  
  unnorm_posterior_densities <- unnorm_likelihoods
  
  norm <- sum(unnorm_posterior_densities) * delta
  posterior_densities <- unnorm_posterior_densities / norm
  
  unif_posterior_means[n] <- sum(xs * posterior_densities) * delta
  unif_posterior_sds[n] <- sqrt(sum( (xs - unif_posterior_means[n])**2 * posterior_densities) * delta)
} 

outer_shrinkage <- data.frame(likelihood_peaks, unif_posterior_means, unif_posterior_sds)

beta1 <- 1
beta2 <- 10
gamma <- 0.75

prior_densities <- gamma * exp(-0.5 * (xs / beta1)**2) / sqrt(2 * pi * beta1**2) +
  (1 - gamma) * exp(-0.5 * (xs / beta2)**2) / sqrt(2 * pi * beta2**2) 

posterior_means <- rep(0, N)
posterior_sds <- rep(0, N)

for (n in 1:N) {
  unnorm_likelihoods <- exp(- 0.5 * ( (xs - likelihood_peaks[n]) / sigma)**2)
  
  unnorm_posterior_densities <- unnorm_likelihoods * prior_densities
  
  norm <- sum(unnorm_posterior_densities) * delta
  posterior_densities <- unnorm_posterior_densities / norm
  
  posterior_means[n] <- sum(xs * posterior_densities) * delta
  posterior_sds[n] <- sqrt(sum( (xs - posterior_means[n])**2 * posterior_densities) * delta)
}

plot(likelihood_peaks, posterior_means / unif_posterior_means, 
     type="l", lwd=2, col=c_dark, xlim=c(5 * delta, 200), ylim=c(0, 1.1))
lines(c(0, 80), c(1, 1), lwd=2, col=c_light)
lines(c(beta1, beta1), c(0, 1.1), lwd=2, col=c_light)
lines(c(beta2, beta2), c(0, 1.1), lwd=2, col=c_light)

plot(likelihood_peaks, posterior_sds / unif_posterior_sds, 
     type="l", lwd=2, col=c_dark, xlim=c(0, 200), ylim=c(0, 1.35))
lines(c(0, 8), c(1, 1), lwd=2, col=c_light)
lines(c(beta, beta), c(0, 1.1), lwd=2, col=c_light)
lines(c(beta2, beta2), c(0, 1.1), lwd=2, col=c_light)

outer_shrinkage['normal_mixture_means'] <- posterior_means
outer_shrinkage['normal_mixture_sds'] <- posterior_sds

write.csv(outer_shrinkage, "outer_shrinkage.csv", row.names=FALSE, quote=FALSE)

############################################################
#
# Outer Extent Density Regularization
#
############################################################

delta <- 0.025
xs <- seq(-5, 125, delta)

density_comps <- data.frame(xs)

# Degenerate Likelihood Function
sigma <- 30
y <- 7
alpha <- 3

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2) / (1 + exp(-alpha * (xs - y)))
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

density_comps['likelihoods'] = likelihoods

############################################################
# Normal Mixture
############################################################

beta1 <- 1
beta2 <- 10
gamma <- 0.75

prior_densities <- gamma * exp(-0.5 * (xs / beta1)**2) / sqrt(2 * pi * beta1**2) +
  (1 - gamma) * exp(-0.5 * (xs / beta2)**2) / sqrt(2 * pi * beta2**2) 

density_comps['normal_mixture_prior'] = prior_densities

unnorm_posterior_densities <- unnorm_likelihoods * prior_densities
norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-5, 125), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, likelihoods, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

density_comps['normal_mixture_posterior'] = posterior_densities

############################################################
# Laplace
############################################################

beta <- 1

prior_densities <- exp( - abs(xs) / beta) / (2 * beta) 

density_comps['laplace_prior'] = prior_densities

unnorm_posterior_densities <- unnorm_likelihoods * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-5, 125), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, likelihoods, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

density_comps['laplace_posterior'] = posterior_densities

############################################################
# Cauchy
############################################################

beta <- 1

prior_densities <- (beta**2 / (xs**2 + beta**2)) / (pi * beta) 

density_comps['cauchy_prior'] = prior_densities

unnorm_posterior_densities <- unnorm_likelihoods * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-5, 125), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, likelihoods, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

density_comps['cauchy_posterior'] = posterior_densities

############################################################
# Save Output
############################################################

write.csv(density_comps, "outer_density_comps.csv", row.names=FALSE, quote=FALSE)


############################################################
#
# Marginal Population Model Comparison
#
############################################################

delta <- 0.025
xs <- seq(-100, 100, delta)

density_comps <- data.frame(xs)

# Normal
beta <- 1
prior_densities <- exp(- 0.5 * (xs / beta)**2) / sqrt(2 * pi * beta**2)

density_comps['normal_prior'] = prior_densities

# Normal Mixture
beta1 <- 1
beta2 <- 10
gamma <- 0.75

prior_densities <- gamma * exp(-0.5 * (xs / beta1)**2) / sqrt(2 * pi * beta1**2) +
  (1 - gamma) * exp(-0.5 * (xs / beta2)**2) / sqrt(2 * pi * beta2**2) 

density_comps['normal_mixture_prior'] = prior_densities

# Laplace
beta <- 1

prior_densities <- exp( - abs(xs) / beta) / (2 * beta) 

density_comps['laplace_prior'] = prior_densities

# Cauchy
beta <- 1

prior_densities <- (beta**2 / (xs**2 + beta**2)) / (pi * beta) 

density_comps['cauchy_prior'] = prior_densities

# Horseshoe
beta <- 1

thetas <- sapply(1:1000000, function(n) rnorm(1, 0, abs(rcauchy(1, 0, 1)) * tau))

density_comps['horseshoe_prior'] <- density(thetas, n=length(xs), from=-100, to=100)$y

############################################################
# Save Output
############################################################

plot(xs, density_comps$normal_prior, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-25, 25), xlab='theta', ylim=c(0, 1), ylab="")

lines(xs, density_comps$normal_mixture_prior, lwd=2, col=c_mid)
lines(xs, density_comps$laplace_prior, lwd=2, col=c_mid_highlight)
lines(xs, density_comps$cauchy_prior, lwd=2, col=c_dark)
lines(xs, density_comps$horseshoe_prior, lwd=2, col=c_dark_highlight)

plot(xs, log(density_comps$normal_prior / density_comps$cauchy_prior), type="l", lwd=2, col=c_light_highlight,
     xlim=c(-25, 25), xlab='theta', ylim=c(-10, 100), ylab="")

lines(xs, log(density_comps$normal_mixture_prior / density_comps$cauchy_prior), lwd=2, col=c_mid)
lines(xs, log(density_comps$laplace_prior / density_comps$cauchy_prior), lwd=2, col=c_mid_highlight)
lines(xs, log(density_comps$cauchy_prior / density_comps$cauchy_prior), lwd=2, col=c_dark)
lines(xs, log(density_comps$horseshoe_prior / density_comps$cauchy_prior), lwd=2, col=c_dark_highlight)

plot(xs, density_comps$horseshoe_prior / density_comps$cauchy_prior, type="l", lwd=2, col=c_dark,
     xlim=c(-50, 50), xlab='theta', ylim=c(0, 3), ylab="")
lines(c(-50, 50), c(1, 1), lwd=2, col=c_light)

write.csv(density_comps, "prior_density_comps.csv", row.names=FALSE, quote=FALSE)


