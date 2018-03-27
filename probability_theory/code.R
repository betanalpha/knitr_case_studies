#setwd('')

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

# Poisson PMF
p <- hist(0, breaks=0:21-0.5, plot=FALSE)
p$counts <- dpois(0:20, 5)

par(mar = c(8, 6, 0, 0.5)) #b-l-t-r
plot(p, main="", col="white", border=c_dark_highlight,
     xlab="x", yaxt='n', ylim=c(0, 0.2), ylab="Probability Mass",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

# Poisson Probability
p_sum <- hist(5, breaks=5:11-0.5, plot=FALSE)
p_sum$counts <- dpois(5:10, 5)

plot(p_sum, col=c_dark, border=c_dark_highlight, add=T)

par(xpd=NA)
arrows(4.5, -0.05, 10.5, -0.05, col = "black", lty = 1, lwd = 2, code=3, angle=90, length=0.1)
text(x=7.5, y=-0.06, labels="A", cex=1.5)
text(x=7, y=0.025, labels="P[A]", cex=1.5, col="white")

# Gaussian PDF
x <- seq(-6, 6, 0.001)

par(mar = c(8, 6, 0, 0.5)) #b-l-t-r
plot(x, dnorm(x, 0, 1), type="l", col=c_dark_highlight, lwd=2,
     xlab="x", ylab="Probability Density",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, yaxt='n')

# Gaussian Probabilty
x_int <- seq(-2, 1, 0.001)
x <- c(x_int, 1, -2)
y <- c(dnorm(x_int, 0, 1), 0, 0)

polygon(x, y, col=c_dark, border=NA)

par(xpd=NA)
arrows(-2, -0.1, 1, -0.1, col = "black", lty = 1, lwd = 2, code=3, angle=90, length=0.1)
text(x=-0.5, y=-0.11, labels="A", cex=1.5)
text(x=-0.5, y=0.1, labels="P[A]", cex=1.5, col="white")

# Poisson Running Mean
N <- 1000

set.seed(13703599)
S <- rpois(N, 5)

emp_exp <- sapply(1:N, function(n) mean(S[1:n]))
emp_var <- sapply(1:N, function(n) var(S[1:n]))
mc_se <- sapply(1:N, function(n) sqrt(var(S[1:n]) / n))

plot(1:N, emp_exp, type="l", col=c_dark_highlight, lwd=2,
ylim=c(2, 6), xlab="N", ylab="Empirical Mean",
cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

abline(h=5, col="grey", lty="dashed", lwd=2)
polygon(c(1:N, N:1), c(emp_exp + mc_se, rev(emp_exp - mc_se)), col=c_mid, border=NA)
lines(1:N, emp_exp, type="l", col=c_dark_highlight, lwd=2)

# Poisson Interval Probability
S_I <- 5 <= S & S <= 10

emp_exp <- sapply(1:N, function(n) mean(S_I[1:n]))
emp_var <- sapply(1:N, function(n) var(S_I[1:n]))
mc_se <- sapply(1:N, function(n) sqrt(var(S_I[1:n]) / n))

plot(1:N, emp_exp, type="l", col=c_dark_highlight, lwd=2,
     ylim=c(2, 6), xlab="N", ylab="Empirical Mean",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

abline(h=ppois(K2, 5) - ppois(K1 - 1, 5), col="grey", lty="dashed", lwd=2)
polygon(c(1:N, N:1), c(emp_exp + mc_se, rev(emp_exp - mc_se)), col=c_mid, border=NA)
lines(1:N, emp_exp, type="l", col=c_dark_highlight, lwd=2)


# Gaussian Running Mean (R)
N <- 1000

S <- rnorm(N, 0, 1)

emp_exp <- sapply(1:N, function(n) mean(S[1:n]))
emp_var <- sapply(1:N, function(n) var(S[1:n]))
mc_se <- sapply(1:N, function(n) sqrt(var(S[1:n]) / n))

plot(1:N, emp_exp, type="l", col=c_dark_highlight, lwd=2,
     ylim=c(-1, 1), xlab="N", ylab="Empirical Mean",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

abline(h=5, col="grey", lty="dashed", lwd=2)
polygon(c(1:N, N:1), c(emp_exp + mc_se, rev(emp_exp - mc_se)), col=c_mid, border=NA)
lines(1:N, emp_exp, type="l", col=c_dark_highlight, lwd=2)

# Gaussian Running Mean (Stan Exact)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source('stan_utility.R')

fit <- stan(file='generate_gaussian.stan', iter=N, warmup=0
            chains=1, seed=194838, algorithm="Fixed_param")

S <- extract(fit)$x[1,]

emp_exp <- sapply(1:N, function(n) mean(S[1:n]))
emp_var <- sapply(1:N, function(n) var(S[1:n]))
mc_se <- sapply(1:N, function(n) sqrt(var(S[1:n]) / n))

plot(1:N, emp_exp, type="l", col=c_dark_highlight, lwd=2,
ylim=c(-1, 1), xlab="N", ylab="Empirical Mean",
cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

abline(h=5, col="grey", lty="dashed", lwd=2)
polygon(c(1:N, N:1), c(emp_exp + mc_se, rev(emp_exp - mc_se)), col=c_mid, border=NA)
lines(1:N, emp_exp, type="l", col=c_dark_highlight, lwd=2)

# Gaussian Running Mean (Stan MCMC)
fit <- stan(file='fit_gaussian.stan', seed=4938483)

source('stan_utility.R')
check_all_diagnostics(fit)

fit_summary <- summary(fit, probs = c(0.5))$summary

S <- extract(fit)$x[1,]

print(fit)
