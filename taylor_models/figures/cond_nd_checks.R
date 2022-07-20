set.seed(2958481)

##################################################
# Prior Pushforward
##################################################

# Prior samples
S <- 1000
prior_beta0 <- rnorm(S, 0, 0.5 / 2.32)
prior_beta1 <- rnorm(S, 0, (0.5 / 2) / 2.32)
prior_beta2 <- rnorm(S, 0, (0.5 / 4) / 2.32)

# Covariate samples and binning
N <- 30
x_sample <- rnorm(N, 0, 1.2)

bins <- seq(-3, 3, 1.5)
B <- length(bins) - 1
bin_idx <- sapply(1:B, function(b) which(bins[b] < x_sample & x_sample < bins[b + 1]))

# Info for a single sample 
idx <- 1
prior_beta0[idx]
prior_beta1[idx]
prior_beta2[idx]

mu_sample <- prior_beta0[idx] + 
             prior_beta1[idx] * x_sample + 
             prior_beta2[idx] * x_sample**2

plot(x_sample, mu_sample, xlim=c(-3, 3), ylim=c(-2, 2))

cat(sprintf("%.3f/%.3f,", x_sample, mu_sample), "\n")

bin_medians <- sapply(1:B, function(b) median(mu_sample[bin_idx[[b]]]))

cat(sprintf("%.3f/%.3f/%.3f,", head(bins, 4), tail(bins, 4), bin_medians), "\n")

# Info for ensemble
mu_samples <- sapply(1:S, function(idx) prior_beta0[idx] + 
                                        prior_beta1[idx] * x_sample +
                                        prior_beta2[idx] * x_sample**2)

median_samples <- sapply(1:S, function(idx) 
                              sapply(1:B, function(b) 
                                          median(mu_samples[bin_idx[[b]], idx])))

cat(sprintf("%.3f/%.3f/%.3f/%.3f,", 
            median_samples[1, 1:25], median_samples[2, 1:25], 
            median_samples[3, 1:25], median_samples[4, 1:25]), "\n")

probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

quantiles <- sapply(1:B, function(b) quantile(median_samples[b,], probs=probs))

r1 <- c(quantiles[1,], rev(quantiles[9,]))
r2 <- c(quantiles[2,], rev(quantiles[8,]))
r3 <- c(quantiles[3,], rev(quantiles[7,]))
r4 <- c(quantiles[4,], rev(quantiles[6,]))
rs <- rbind(r1, r2, r3, r4)

cat(sprintf("%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f,", 
            rs[,1], rs[,2], rs[,3], rs[,4], 
            rs[,5], rs[,6], rs[,7], rs[,8]), "\n")

quantiles[5,]

plot(1, type="n", main="",
     xlim=c(-3, 3), xlab="x", ylim=c(-2, 2), ylab="")
points(xs, r1, col=c_light, pch=16)
points(xs, r2, col=c_light_highlight, pch=16)
points(xs, r3, col=c_mid, pch=16)
points(xs, r4, col=c_mid_highlight, pch=16)

##################################################
# Posterior Predictive Pushforward
##################################################

# Observed data
idx <- 1
true_beta0 <- prior_beta0[idx]
true_beta1 <- prior_beta1[idx]
true_beta2 <- prior_beta2[idx]

obs <- true_beta0 + 
       true_beta1 * x_sample + 
       true_beta2 * x_sample**2 + 
       rnorm(N, 0, 0.5)

plot(x_sample, obs, xlim=c(-3, 3), ylim=c(-2, 2))

cat(sprintf("%.3f/%.3f,", x_sample, obs), "\n")

obs_bin_medians <- sapply(1:B, function(b) median(obs[bin_idx[[b]]]))

cat(sprintf("%.3f/%.3f/%.3f,", head(bins, 4), tail(bins, 4), obs_bin_medians), "\n")

# Posterior predictions
post_beta0 <- rnorm(S, true_beta0, 0.01)
post_beta1 <- rnorm(S, true_beta1, 0.01)
post_beta2 <- rnorm(S, true_beta2, 0.01)

cat(sprintf("%.3f/%.3f/%.3f,", post_beta0[1:20], post_beta1[1:20], post_beta2[1:20]), "\n")

x_grid <- seq(-3, 3, 0.05)
G <- length(x_grid)

grid_pred_samples <- sapply(1:S, function(idx) post_beta0[idx] + 
                            post_beta1[idx] * x_grid +
                            post_beta2[idx] * x_grid**2 +
                            rnorm(G, 0, 0.5))

quantiles <- sapply(1:G, function(g) quantile(grid_pred_samples[g,], probs=probs))


cat(sprintf("%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f,", 
            x_grid, quantiles[1,], quantiles[2,], quantiles[3,], 
            quantiles[4,], quantiles[5,], quantiles[6,],
            quantiles[7,], quantiles[8,], quantiles[9,] ), "\n")


pred_samples <- sapply(1:S, function(idx) post_beta0[idx] + 
                                          post_beta1[idx] * x_sample +
                                          post_beta2[idx] * x_sample**2 +
                                          rnorm(N, 0, 0.5))



bin_medians <- sapply(1:S, function(idx) 
  sapply(1:B, function(b) 
    median(pred_samples[bin_idx[[b]], idx])))


cat(sprintf("%.3f/%.3f,", x_sample, pred_samples[,1]), "\n")
cat(sprintf("%.3f/%.3f/%.3f,", head(bins, 4), tail(bins, 4), bin_medians[,1]), "\n")

cat(sprintf("%.3f/%.3f/%.3f/%.3f,", 
            bin_medians[1, 1:25], bin_medians[2, 1:25], 
            bin_medians[3, 1:25], bin_medians[4, 1:25]), "\n")

quantiles <- sapply(1:B, function(b) quantile(bin_medians[b,], probs=probs))

r1 <- c(quantiles[1,], rev(quantiles[9,]))
r2 <- c(quantiles[2,], rev(quantiles[8,]))
r3 <- c(quantiles[3,], rev(quantiles[7,]))
r4 <- c(quantiles[4,], rev(quantiles[6,]))
rs <- rbind(r1, r2, r3, r4)

cat(sprintf("%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f,", 
            rs[,1], rs[,2], rs[,3], rs[,4], 
            rs[,5], rs[,6], rs[,7], rs[,8]), "\n")

quantiles[5,]


bin_centers <- sapply(1:B, function(b) 0.5 * (bins[b] + bins[b + 1]))
xs <- c(bin_centers, rev(bin_centers))

plot(1, type="n", main="",
     xlim=c(-3, 3), xlab="x", ylim=c(-2, 2), ylab="")
points(xs, r1, col=c_light, pch=16)
points(xs, r2, col=c_light_highlight, pch=16)
points(xs, r3, col=c_mid, pch=16)
points(xs, r4, col=c_mid_highlight, pch=16)
points(bin_centers, obs_bin_medians, col="black", pch=16)

range(pred_samples)

# Marginal
marginal_bins <- seq(-3, 3, 0.5)
mB <- length(marginal_bins) - 1

idx <- rep(1:mB, each=2)
xs <- sapply(1:length(idx), function(b) if(b %% 2 == 0) marginal_bins[idx[b] + 1] else marginal_bins[idx[b]] )

obs_counts <- hist(obs, breaks=marginal_bins)$counts
pad_obs_counts <- obs_counts[idx]

cat(sprintf("%.3f/%.3f,", xs, pad_obs_counts), "\n")
    
post_pred_counts <- sapply(1:S, function(s) 
                                hist(pred_samples[,s], 
                                     breaks=marginal_bins, plot=FALSE)$counts)

cred <- sapply(1:mB, function(b) quantile(post_pred_counts[b,], probs=probs))

cat(sprintf("%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f/%.3f,", 
            cred[1,], cred[2,], cred[3,], cred[4,], cred[5,],
            cred[6,], cred[7,], cred[8,], cred[9,], 
            head(marginal_bins, mB), tail(marginal_bins, mB)), "\n")


# Failure One
idx <- 1
true_beta0 <- prior_beta0[idx]
true_beta1 <- prior_beta1[idx]
true_beta2 <- prior_beta2[idx]

obs <- true_beta0 + 
       true_beta1 * x_sample + 
       true_beta2 * x_sample**2 + 
       0.15 * x_sample**3 +
       rnorm(N, 0, 0.5)

plot(x_sample, obs, xlim=c(-3, 3), ylim=c(-2, 2))

cat(sprintf("%.3f/%.3f,", x_sample, obs), "\n")

marginal_bins <- seq(-3, 3, 0.5)
mB <- length(marginal_bins) - 1

idx <- rep(1:mB, each=2)
xs <- sapply(1:length(idx), function(b) if(b %% 2 == 0) marginal_bins[idx[b] + 1] else marginal_bins[idx[b]] )

obs_counts <- hist(obs, breaks=marginal_bins)$counts
pad_obs_counts <- obs_counts[idx]

cat(sprintf("%.3f/%.3f,", xs, pad_obs_counts), "\n")

# Failure Two
idx <- 1
true_beta0 <- prior_beta0[idx]
true_beta1 <- prior_beta1[idx]
true_beta2 <- prior_beta2[idx]

obs <- true_beta0 + 
  true_beta1 * x_sample + 
  true_beta2 * x_sample**2 + 
  0.5 * rt(N, 3)

plot(x_sample, obs, xlim=c(-3, 3), ylim=c(-2, 2))

cat(sprintf("%.3f/%.3f,", x_sample, obs), "\n")

marginal_bins <- seq(-3, 3, 0.5)
mB <- length(marginal_bins) - 1

idx <- rep(1:mB, each=2)
xs <- sapply(1:length(idx), function(b) if(b %% 2 == 0) marginal_bins[idx[b] + 1] else marginal_bins[idx[b]] )

obs_counts <- hist(obs, breaks=marginal_bins)$counts
pad_obs_counts <- obs_counts[idx]

cat(sprintf("%.3f/%.3f,", xs, pad_obs_counts), "\n")

# Failure Two
idx <- 1
true_beta0 <- prior_beta0[idx]
true_beta1 <- prior_beta1[idx]
true_beta2 <- prior_beta2[idx]

obs <- true_beta0 + 
  true_beta1 * x_sample + 
  true_beta2 * x_sample**2 + 
  0.5 * rt(N, 3)

plot(x_sample, obs, xlim=c(-3, 3), ylim=c(-2, 2))

cat(sprintf("%.3f/%.3f,", x_sample, obs), "\n")

marginal_bins <- seq(-3, 3, 0.5)
mB <- length(marginal_bins) - 1

idx <- rep(1:mB, each=2)
xs <- sapply(1:length(idx), function(b) if(b %% 2 == 0) marginal_bins[idx[b] + 1] else marginal_bins[idx[b]] )

obs_counts <- hist(obs, breaks=marginal_bins)$counts
pad_obs_counts <- obs_counts[idx]

cat(sprintf("%.3f/%.3f,", xs, pad_obs_counts), "\n")



# Correlated Covarites
N <- 50
x1 <- rnorm(N, 0, 1.25)
x2 <- rnorm(N, x1**2, 0.5)

plot(x1, x2)

cat(sprintf("%.3f/%.3f,", x1, x2), "\n")

hist(x1, breaks=bins)$counts

idx <- 1
true_beta0 <- prior_beta0[idx]
true_beta1 <- prior_beta1[idx]
true_beta2 <- prior_beta2[idx]

obs <- true_beta0 + 
  true_beta1 * x1 + 
  true_beta2 * x1**2 + 
  rnorm(N, 0, 0.5)

plot(x1, obs, xlim=c(-3, 3), ylim=c(-2, 2))

cat(sprintf("%.3f/%.3f,", x_sample, obs), "\n")

bins <- seq(-3, 3, 0.5)
B <- length(bins) - 1
bin_idx <- sapply(1:B, function(b) which(bins[b] < x_sample & x_sample < bins[b + 1]))

bin_biases <- sapply(1:B, function(b) (0.5 *(bins[b] + bins[b + 1]))**2)

obs_bin_medians <- sapply(1:B, function(b) median(obs[bin_idx[[b]]]) )

cat(sprintf("%.3f/%.3f/%.3f,", head(bins, B), tail(bins, B), 
                               obs_bin_medians + bin_biases), "\n")

post_beta0 <- rnorm(S, true_beta0, 0.01)
post_beta1 <- rnorm(S, true_beta1, 0.01)
post_beta2 <- rnorm(S, true_beta2, 0.01)

pred_samples <- sapply(1:S, function(idx) post_beta0[idx] + 
                         post_beta1[idx] * x1 +
                         post_beta2[idx] * x2**2 +
                         rnorm(N, 0, 0.5))

bin_medians <- sapply(1:S, function(idx) 
  sapply(1:B, function(b) 
    median(pred_samples[bin_idx[[b]], idx])))

bin_medians <- bin_medians + bin_biases

cat(sprintf("%.3f/%.3f,", x_sample, pred_samples[,1]), "\n")
cat(sprintf("%.3f/%.3f/%.3f,", head(bins, 4), tail(bins, 4), bin_medians[,1]), "\n")

cat(sprintf("%.3f/%.3f/%.3f/%.3f,", 
            bin_medians[1, 1:25], bin_medians[2, 1:25], 
            bin_medians[3, 1:25], bin_medians[4, 1:25]), "\n")

quantiles <- sapply(1:B, function(b) quantile(bin_medians[b,], probs=probs, na.rm=TRUE))

cat(sprintf("%.3f/%.3f/%.3f/%.3f,", 
            head(bins, B), tail(bins, B), quantiles[1,], quantiles[9,]), "\n")

cat(sprintf("%.3f/%.3f/%.3f/%.3f,", 
            head(bins, B), tail(bins, B), quantiles[2,], quantiles[8,]), "\n")

cat(sprintf("%.3f/%.3f/%.3f/%.3f,", 
            head(bins, B), tail(bins, B), quantiles[3,], quantiles[7,]), "\n")

cat(sprintf("%.3f/%.3f/%.3f/%.3f,", 
            head(bins, B), tail(bins, B), quantiles[4,], quantiles[6,]), "\n")

cat(sprintf("%.3f/%.3f/%.3f,", 
            head(bins, B), tail(bins, B), quantiles[5,]), "\n")


cat(sprintf("%.3f/%.3f/%.3f,", 
            head(bins, B), tail(bins, B), obs_bin_medians + bin_biases - quantiles[5,]), "\n")

cat(sprintf("%.3f/%.3f/%.3f/%.3f,", 
            head(bins, B), tail(bins, B), quantiles[1,] - quantiles[5,], quantiles[9,] - quantiles[5,]), "\n")

cat(sprintf("%.3f/%.3f/%.3f/%.3f,", 
            head(bins, B), tail(bins, B), quantiles[2,] - quantiles[5,], quantiles[8,] - quantiles[5,]), "\n")

cat(sprintf("%.3f/%.3f/%.3f/%.3f,", 
            head(bins, B), tail(bins, B), quantiles[3,] - quantiles[5,], quantiles[7,] - quantiles[5,]), "\n")

cat(sprintf("%.3f/%.3f/%.3f/%.3f,", 
            head(bins, B), tail(bins, B), quantiles[4,] - quantiles[5,], quantiles[6,] - quantiles[5,]), "\n")

cat(sprintf("%.3f/%.3f/%.3f,", 
            head(bins, B), tail(bins, B), 0), "\n")