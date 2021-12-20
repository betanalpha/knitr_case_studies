library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

data <- read_rdump("poly.data.R")

set.seed(4838583)

# Prior
fit <- stan(file='poly_prior.stan', data=data, 
            seed=5838299, refresh=0)

samples <- extract(fit)

library(colormap)

I <- 4000

plot_idx <- seq(1, I, 50)
J <- length(plot_idx)

nom_colors <- c("#DCBCBC", "#C79999", "#B97C7C",
                "#A25050", "#8F2727", "#7C0000")
line_colors <- colormap(colormap=nom_colors, nshades=2 * J)

plot(1, type="n", xlab="x", ylab="y", 
     main="",
     xlim=c(-11, 11), ylim=c(-75, 75))
for (j in 1:J)
  lines(data$x_predict, samples$f[plot_idx[j],], col=line_colors[0.25 * J + j], lwd=2)

write.csv(data.frame(samples$f), file="output_prior.clean.csv", row.names=F)

probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:data$N_predict, 
               function(n) quantile(samples$f[,n], probs=probs))

plot(1, type="n", main="",
     xlab="x",  ylab="y", xlim=c(-11, 11), ylim=c(-75, 75))
polygon(c(data$x_predict, rev(data$x_predict)), c(cred[1,], rev(cred[9,])),
        col = c_light, border = NA)
polygon(c(data$x_predict, rev(data$x_predict)), c(cred[2,], rev(cred[8,])),
        col = c_light_highlight, border = NA)
polygon(c(data$x_predict, rev(data$x_predict)), c(cred[3,], rev(cred[7,])),
        col = c_mid, border = NA)
polygon(c(data$x_predict, rev(data$x_predict)), c(cred[4,], rev(cred[6,])),
        col = c_mid_highlight, border = NA)
lines(data$x_predict, cred[5,], col=c_dark, lwd=2)

df <- data.frame(t(cred))
df <- cbind(data$x_predict, df)

write.csv(df, file="cred_prior.csv", row.names=F)

# Posterior
fit <- stan(file='poly.stan', data=data, 
            seed=5838299, refresh=0)

samples <- extract(fit)

plot(1, type="n", xlab="x", ylab="y", 
     main="",
     xlim=c(-11, 11), ylim=c(-75, 75))
for (j in 1:J)
  lines(data$x_predict, samples$f[plot_idx[j],], col=line_colors[0.25 * J + j], lwd=2)

write.csv(data.frame(samples$f), file="output.clean.csv", row.names=F)

cred <- sapply(1:data$N_predict, 
               function(n) quantile(samples$f[,n], probs=probs))

plot(1, type="n", main="",
     xlab="x",  ylab="y", xlim=c(-11, 11), ylim=c(-75, 75))
polygon(c(data$x_predict, rev(data$x_predict)), c(cred[1,], rev(cred[9,])),
        col = c_light, border = NA)
polygon(c(data$x_predict, rev(data$x_predict)), c(cred[2,], rev(cred[8,])),
        col = c_light_highlight, border = NA)
polygon(c(data$x_predict, rev(data$x_predict)), c(cred[3,], rev(cred[7,])),
        col = c_mid, border = NA)
polygon(c(data$x_predict, rev(data$x_predict)), c(cred[4,], rev(cred[6,])),
        col = c_mid_highlight, border = NA)
lines(data$x_predict, cred[5,], col=c_dark, lwd=2)

df <- data.frame(t(cred))
df <- cbind(data$x_predict, df)

write.csv(df, file="cred.csv", row.names=F)

