alpha <- -0.5
beta <- 0.75

N <- 5
xs <- rnorm(N, 0, 1)
ys <- rnorm(N, alpha + beta * xs, 0.25)

plot(xs, ys, pch=18, xlim=c(-3, 3), ylim=c(-2, 2))

X <- cbind(rep(1, N), xs)
Z <- solve(t(X) %*% X)

mu <- ys %*% X %*% Z
Sigma <- (0.25)**2 * Z
L <- chol(Sigma)

sigma_x <- sqrt(Sigma[1,1])
sigma_y <- sqrt(Sigma[2,2])
rho <- Sigma[1, 2] / (sigma_x * sigma_y)

sigma_x; sigma_y; rho


plot(xs, ys, pch=18, xlim=c(-3, 3), ylim=c(-2, 2))

xps <- seq(-3, 3, 0.01)


gammas <- matrix(NA, nrow=2, ncol=10)

for (n in 1:10) {
  zs <- rnorm(2, 0, 1)
  gamma <- t(mu) + t(L) %*% zs
  lines(xps, gamma[1] + gamma[2] * xps, col="red")
  gammas[,n] <- gamma
}

cat(sprintf("%.3f/%.3f,", xs, ys), "\n")
cat(sprintf("%.3f/%.3f,", gammas[1,], gammas[2,]), "\n")

x <- xs[3]
y <- ys[3]

plot(x, y, pch=18, xlim=c(-3, 3), ylim=c(-2, 2))

gammas <- matrix(NA, nrow=2, ncol=10)

for (n in 1:10) {
  alpha_nom <- rnorm(1, 0, 1)
  beta_nom <- ((y - alpha_nom) / x)
  
  delta <- rnorm(1, 0, 0.5)
  alpha <- alpha_nom + delta
  beta <- beta_nom + (x * alpha_nom) * delta
  
  lines(xps, alpha + beta * xps, col="red")
  gammas[,n] <- c(alpha, beta)
}

cat(sprintf("%.3f/%.3f,", gammas[1,], gammas[2,]), "\n")


