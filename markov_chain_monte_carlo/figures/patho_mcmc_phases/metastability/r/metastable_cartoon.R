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

############################################################
#
# How does the Random Walk Metropolis algorithm perform
# on a target distribution with a two-dimensional Gaussian 
# density function?
#
############################################################

D <- 2

target_lpdf1 <- function(x) {
  xc = x[1] + 5
  yc = x[2] - 2
  r = sqrt(xc**2 + yc**2)
  a = 3
  b = 5
  e2 = 1 - (b / a)**2
  cos2t = yc**2 / (xc**2 + yc**2)
  r0 = b / sqrt(1 - e2 * cos2t)
  sigma = 0.2
  (-0.5 * ((r - r0) / sigma)**2 - 0.5 * log(6.283185307179586) - log(sigma))
}

target_lpdf2 <- function(x) {
  xc = x[1] - 5
  yc = x[2] + 4.5
  r = sqrt(xc**2 + yc**2)
  a = 3
  b = 5
  e2 = 1 - (b / a)**2
  cos2t = yc**2 / (xc**2 + yc**2)
  r0 = b / sqrt(1 - e2 * cos2t)
  sigma = 0.2
  (-0.5 * ((r - r0) / sigma)**2 - 0.5 * log(6.283185307179586) - log(sigma))
}

target_lpdf <- function(x) {
  lpdf1 <- target_lpdf1(x)
  lpdf2 <- target_lpdf2(x)
  if (lpdf1 > lpdf2) {
    lpdf <- lpdf1 + log(1 + exp(lpdf2 - lpdf1))
  } else {
    lpdf <- lpdf2 + log(1 + exp(lpdf1 - lpdf2))
  }
  (lpdf)
}

N <- 200
q.1 <- seq(-12, 12, 24 / N)
q.2 <- seq(-9.5, 6.5, 16 / N)
densities <- matrix(data = 0, nrow = N + 1, ncol = N + 1)

for (n in 1:N) {
  for (m in 1:N) {
    q <- c(q.1[n], q.2[m])
    densities[n, m] <- exp(target_lpdf(q)) 
  }
}

# To ensure accurate results let's generate pretty large samples
n_init <- 175
n_total <- 200
n_final <- n_total - n_init - 2
n_thin <- 200

# Tune proposal density
sigma <- 0.5

set.seed(8675313)
mcmc_samples <- matrix(data = 0, nrow = n_total, ncol = D)
  
# Randomly seed the initial state
mcmc_samples[1, 1:D] <- c(-8, 5)
  
for (n in 1:n_init) {
  x0 <- mcmc_samples[n, 1:D]
  for (t in 1:n_thin) {
    xp <- rnorm(D, x0, sigma)

    accept_prob <- min(1, exp(target_lpdf(xp) - target_lpdf(x0)))
  
    u = runif(1, 0, 1)
    if (accept_prob > u)
      x0 <- xp
  }
  mcmc_samples[n + 1, 1:D] <- x0
}

mcmc_samples[n_init + 2, 1:D] <- c(2.5, -1.75)

for (n in 1:n_final) {
  x0 <- mcmc_samples[n + n_init + 1, 1:D]
  for (t in 1:n_thin) {
    xp <- rnorm(D, x0, sigma)
    
    accept_prob <- min(1, exp(target_lpdf(xp) - target_lpdf(x0)))
    
    u = runif(1, 0, 1)
    if (accept_prob > u)
      x0 <- xp
  }
  mcmc_samples[n + n_init + 2, 1:D] <- x0
}

contour(q.1, q.2, densities, nlevels=25, main="Gaussian Target Density",
        xlab="q.1", xlim=c(-12, 12), 
        ylab="q.2", ylim=c(-9.5, 6.5), drawlabels = FALSE, col = c_dark, lwd = 2)

points(mcmc_samples[,1], mcmc_samples[,2])
points(mcmc_samples[175, 1], mcmc_samples[175, 2], col="green")

distances <- sqrt((mcmc_samples[,1] + 1.5)**2 + (mcmc_samples[,2] + 0.25)**2)
plot(150:180, distances[150:180])

plot(1:n_total, mcmc_samples[,1], xlim=c(170,180))
plot(1:n_total, mcmc_samples[,2], xlim=c(170,180))

plot(1:n_total, cumsum(mcmc_samples[,2]) / 1:n_total)

for (n in 1:n_total) {
  print(paste("mcpoint{", mcmc_samples[n,1], "}{", mcmc_samples[n,2], "}", sep=""))
}
     
for (n in 1:n_total) {
  print(c(mcmc_samples[n,1], mcmc_samples[n,2]))
}