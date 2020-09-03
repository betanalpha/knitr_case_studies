library(colormap)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

nom_colors <- c("#DCBCBC", "#C79999", "#B97C7C", "#A25050", "#8F2727", "#7C0000")

c_light_teal <- c("#6B8E8E")
c_mid_teal <- c("#487575")
c_dark_teal <- c("#1D4F4F")


# Plot Gaussian process realizations
plot_gp_prior_realizations <- function(fit, xs, title) {
  samples <- extract(fit)
  I <- length(samples$f[,1])
  
  plot_idx <- seq(1, I, 80)
  N <- length(plot_idx)
  line_colors <- colormap(colormap=nom_colors, nshades=N)

  plot(1, type="n", xlab="x", ylab="f", main=title,
       xlim=c(-11, 11), ylim=c(-12, 12))
  for (n in 1:N)
    lines(xs, samples$f[plot_idx[n],], col=line_colors[n], lwd=2)
}

plot_gp_post_realizations <- function(fit, data, true, title) {
  samples <- extract(fit)
  I <- length(samples$f_predict[,1])

  plot_idx <- seq(1, I, 80)
  N <- length(plot_idx)
  line_colors <- colormap(colormap=nom_colors, nshades=N)

  plot(1, type="n", xlab="x", ylab="f", main=title,
       xlim=c(-11, 11), ylim=c(-12, 12))
  for (n in 1:N)
    lines(data$x_predict, samples$f_predict[plot_idx[n],], col=line_colors[n], lwd=2)
    
  lines(true$x, true$f, lwd=4, xlab="x", ylab="f", col="white")
  lines(true$x, true$f, lwd=2, xlab="x", ylab="f", col="black")
}

# Plot low sigma Gaussian process realizations
plot_low_sigma_gp_post_realizations <- function(fit, data, true, title) {
  samples <- extract(fit)
  I <- length(samples$f_predict[,1])

  plot(1, type="n", xlab="x", ylab="f", main=title,
       xlim=c(-11, 11), ylim=c(-10, 10))

  for (i in seq(1, I, 40))
    lines(data$x_predict, samples$f_predict[i,], col=c("#CCCCCC"), lwd=2)

  plot_idx <- which(samples$sigma < 0.5)
  N <- length(plot_idx)
  c_superfine <- c("#8F272710")
  for (n in 1:N)
    lines(data$x_predict, samples$f_predict[plot_idx[n],], col=c_superfine, lwd=2)

  lines(true$x, true$f, lwd=4, xlab="x", ylab="f", col="white")
  lines(true$x, true$f, lwd=2, xlab="x", ylab="f", col="black")

  points(data$x_obs, data$y_obs, col="white", pch=16, cex=1.2)
  points(data$x_obs, data$y_obs, col="black", pch=16, cex=0.8)
}

# Plot Gaussian process quantiles
plot_gp_prior_quantiles <- function(fit, xs, title) {
  samples <- extract(fit)

  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:length(xs),
                 function(n) quantile(samples$f[,n], probs=probs))

  plot(1, type="n", main=title,
       xlab="x", ylab="f", xlim=c(-11, 11), ylim=c(-12, 12))
  polygon(c(xs, rev(xs)), c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(xs, rev(xs)), c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(xs, rev(xs)), c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(xs, rev(xs)), c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(xs, cred[5,], col=c_dark, lwd=2)
}

plot_gp_post_quantiles <- function(fit, data, true, title) {
  samples <- extract(fit)

  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:length(data$x_predict),
                 function(n) quantile(samples$f_predict[,n], probs=probs))

  plot(1, type="n", main=title,
       xlab="x", ylab="f", xlim=c(-11, 11), ylim=c(-10, 10))
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(data$x_predict, cred[5,], col=c_dark, lwd=2)

  lines(true$x, true$f, lwd=4, xlab="x", ylab="f", col="white")
  lines(true$x, true$f, lwd=2, xlab="x", ylab="f", col="black")
}

# Plot Gaussian process predictive quantiles
plot_gp_prior_pred_quantiles <- function(fit, xs, title) {
  samples <- extract(fit)

  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:length(xs),
                 function(n) quantile(samples$y[,n], probs=probs))

  plot(1, type="n", main=title,
       xlab="x", ylab="y", xlim=c(-11, 11), ylim=c(-12, 12))
  polygon(c(xs, rev(xs)), c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(xs, rev(xs)), c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(xs, rev(xs)), c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(xs, rev(xs)), c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(xs, cred[5,], col=c_dark, lwd=2)
}

plot_gp_prior_pred_quantiles_disc <- function(fit, xs, title) {
  samples <- extract(fit)

  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:length(xs),
                 function(n) quantile(samples$y[,n], probs=probs))

  plot(1, type="n", main=title,
       xlab="x", ylab="y", xlim=c(-11, 11), ylim=c(0, 60))
  polygon(c(xs, rev(xs)), c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(xs, rev(xs)), c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(xs, rev(xs)), c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(xs, rev(xs)), c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(xs, cred[5,], col=c_dark, lwd=2)
}

plot_gp_post_pred_quantiles <- function(fit, data, true, title) {
  samples <- extract(fit)

  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:length(data$x_predict),
                 function(n) quantile(samples$y_predict[,n], probs=probs))

  plot(1, type="n", main=title,
       xlab="x", ylab="y", xlim=c(-11, 11), ylim=c(-10, 10))
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(data$x_predict, cred[5,], col=c_dark, lwd=2)

  points(data$x_predict, data$y_predict, col="white", pch=16, cex=0.6)
  points(data$x_predict, data$y_predict, col=c_mid_teal, pch=16, cex=0.4)
  
  lines(true$x, true$f, lwd=4, xlab="x", ylab="f", col="white")
  lines(true$x, true$f, lwd=2, xlab="x", ylab="f", col="black")
  
  points(data$x_obs, data$y_obs, col="white", pch=16, cex=1.2)
  points(data$x_obs, data$y_obs, col="black", pch=16, cex=0.8)
}

plot_gp_post_pred_quantiles_disc <- function(fit, data, true, title) {
  samples <- extract(fit)

  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:length(data$x_predict),
                 function(n) quantile(samples$y_predict[,n], probs=probs))

  plot(1, type="n", main=title,
       xlab="x", ylab="y", xlim=c(-11, 11), ylim=c(0, 8 ))
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(data$x_predict, rev(data$x_predict)), c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(data$x_predict, cred[5,], col=c_dark, lwd=2)

  points(data$x_predict, data$y_predict, col="white", pch=16, cex=0.6)
  points(data$x_predict, data$y_predict, col=c_mid_teal, pch=16, cex=0.4)

  lines(true$x, true$f, lwd=4, xlab="x", ylab="f", col="white")
  lines(true$x, true$f, lwd=2, xlab="x", ylab="f", col="black")
  
  points(data$x_obs, data$y_obs, col="white", pch=16, cex=1.2)
  points(data$x_obs, data$y_obs, col="black", pch=16, cex=0.8)
}
