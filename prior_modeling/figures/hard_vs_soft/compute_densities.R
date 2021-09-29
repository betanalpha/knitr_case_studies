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
xs <- seq(-20, 20, delta)

density_comps <- data.frame(xs)

# Containment
sigma <- 5
y <- 2

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

density_comps['containment_like'] = likelihoods

# Contraction
sigma <- 0.25
y <- -1.5

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

density_comps['contraction_like'] = likelihoods

# Compromise
sigma <- 0.25
y <- 10.25

unnorm_likelihoods <- exp(- 0.5 * ( (xs - y) / sigma)**2)
likelihoods <- unnorm_likelihoods / (sum(unnorm_likelihoods) * delta)

density_comps['compromise_like'] = likelihoods

############################################################
# Normal
############################################################

beta <- 2

prior_densities <- exp(- 0.5 * (xs / beta)**2) / sqrt(2 * pi * beta**2)

density_comps['normal_prior'] = prior_densities

# Containment
unnorm_posterior_densities <- density_comps$containment_like * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

density_comps['normal_post_containment'] = posterior_densities

par(mfrow=c(1, 1))

plot(xs, prior_densities, type="l", lwd=2, col=c_light,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, density_comps$containment_like, lwd=2, col=c_mid)

lines(xs, posterior_densities, lwd=2, col=c_dark)

# Contraction
unnorm_posterior_densities <- density_comps$contraction_like * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

density_comps['normal_post_contraction'] = posterior_densities

par(mfrow=c(1, 1))

plot(xs, prior_densities, type="l", lwd=2, col=c_light,
     xlim=c(-3, 0), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, density_comps$contraction_like, lwd=2, col=c_mid)

lines(xs, posterior_densities, lwd=2, col=c_dark)

# Compromise
unnorm_posterior_densities <- density_comps$compromise_like * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

density_comps['normal_post_compromise'] = posterior_densities

par(mfrow=c(1, 1))

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, density_comps$compromise_like, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

############################################################
# Uniform
############################################################

beta <- 2

prior_densities <- sapply(xs, function(x) if( abs(x) > beta * 2.57) 0 else 1 / (beta * 2 * 2.57))
  
density_comps['unif_prior'] = prior_densities

# Containment
unnorm_posterior_densities <- density_comps$containment_like * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

density_comps['unif_post_containment'] = posterior_densities

par(mfrow=c(1, 1))

plot(xs, prior_densities, type="l", lwd=2, col=c_light,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, density_comps$containment_like, lwd=2, col=c_mid)

lines(xs, posterior_densities, lwd=2, col=c_dark)

# Contraction
unnorm_posterior_densities <- density_comps$contraction_like * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

density_comps['unif_post_contraction'] = posterior_densities

par(mfrow=c(1, 1))

plot(xs, prior_densities, type="l", lwd=2, col=c_light,
     xlim=c(-3, 0), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, density_comps$contraction_like, lwd=2, col=c_mid)

lines(xs, posterior_densities, lwd=2, col=c_dark)

# Compromise
unnorm_posterior_densities <- density_comps$compromise_like * prior_densities

norm <- sum(unnorm_posterior_densities) * delta
posterior_densities <- unnorm_posterior_densities / norm

density_comps['unif_post_compromise'] = posterior_densities

par(mfrow=c(1, 1))

plot(xs, prior_densities, type="l", lwd=2, col=c_light_highlight,
     xlim=c(-15, 15), xlab='theta', ylim=c(0, 0.5), ylab="")

lines(xs, density_comps$compromise_like, lwd=2, col=c_mid_highlight)

lines(xs, posterior_densities, lwd=2, col=c_dark)

############################################################
# Save Output
############################################################

write.csv(density_comps, "density_comps.csv", row.names=FALSE, quote=FALSE)


