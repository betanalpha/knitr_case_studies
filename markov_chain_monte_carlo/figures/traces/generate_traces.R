N <- 1000

rho <- 0.25

innovations <- rnorm(N, 0, 1)
x <- rep(0, N)

x[1] <- innovations[1]
for (n in 2:N) {
  x[n] <- rho * x[n - 1] + innovations[n]
}

write.table(data.frame(1:N, x), file = "low_rho.dat", row.names=FALSE, col.names=FALSE, sep=" ")

rho <- 0.95

innovations <- rnorm(N, 0, 1)
x <- rep(0, N)

x[1] <- innovations[1]
for (n in 2:N) {
  x[n] <- rho * x[n - 1] + innovations[n]
}

write.table(data.frame(1:N, x), file = "high_rho.dat", row.names=FALSE, col.names=FALSE, sep=" ")

