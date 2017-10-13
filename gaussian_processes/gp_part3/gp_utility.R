c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_trans <- c("#DCBCBC80")
c_light_highlight_trans <- c("#C7999980")
c_mid_trans <- c("#B97C7C80")
c_mid_highlight_trans <- c("#A2505080")
c_dark_trans <- c("#8F272780")
c_dark_highlight_trans <- c("#7C000080")

c_light_teal <- c("#6B8E8E")
c_mid_teal <- c("#487575")
c_dark_teal <- c("#1D4F4F")

# Plot Gaussian process realizations
plot_gp_realizations <- function(fit, data, true, title) {
  params <- extract(fit)
  I <- length(params$f_predict[,1])

  c_superfine <- c("#8F272705")
  plot(1, type="n", xlab="x", ylab="y", main=title,
       xlim=c(-10, 10), ylim=c(-10, 10))
  for (i in 1:I)
     lines(data$x_predict, params$f_predict[i,], col=c_superfine)

  points(data$x_predict, data$y_predict, col="white", pch=16, cex=0.6)
  points(data$x_predict, data$y_predict, col=c_mid_teal, pch=16, cex=0.4)
  lines(true$x, true$f, lwd=2, xlab="x", ylab="y")
  points(data$x, data$y, col="white", pch=16, cex=1.2)
  points(data$x, data$y, col="black", pch=16, cex=0.8)
}

# Plot low sigma Gaussian process realizations
plot_low_sigma_gp_realizations <- function(fit, data, true, title) {
  params <- extract(fit)
  I <- length(params$f_predict[,1])

  c_superfine <- c("#8F272705")
  plot(1, type="n", xlab="x", ylab="y", main=title,
       xlim=c(-10, 10), ylim=c(-10, 10))
  for (i in 1:I)
    lines(data$x_predict, params$f_predict[i,], col=c("#CCCCCC"))
  for (i in 1:I)
    if (params$sigma[i] < 0.5)
      lines(data$x_predict, params$f_predict[i,], col=c_superfine)

  points(data$x_predict, data$y_predict, col="white", pch=16, cex=0.6)
  points(data$x_predict, data$y_predict, col=c_mid_teal, pch=16, cex=0.4)
  lines(true$x, true$f, lwd=2, xlab="x", ylab="y")
  points(data$x, data$y, col="white", pch=16, cex=1.2)
  points(data$x, data$y, col="black", pch=16, cex=0.8)
}

# Plot Gaussian process predictive realizations
plot_gp_pred_realizations <- function(fit, data, true, title) {
  params <- extract(fit)
  I <- length(params$y_predict[,1])

  plot(1, type="n", xlab="x", ylab="y", main=title,
       xlim=c(-10, 10), ylim=c(-10, 10))
  for (i in 1:I)
     lines(data$x_predict, params$y_predict[i,], col=c_superfine)

  points(data$x_predict, data$y_predict, col="white", pch=16, cex=0.6)
  points(data$x_predict, data$y_predict, col=c_mid_teal, pch=16, cex=0.4)
  lines(true$x, true$y, lwd=2, xlab="x", ylab="y")
  points(data$x, data$y, col="white", pch=16, cex=1.2)
  points(data$x, data$y, col="black", pch=16, cex=0.8)
}

# Plot Gaussian process quantiles
plot_gp_quantiles <- function(fit, data, true, title) {
  params <- extract(fit)

  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:length(data$x_predict),
                 function(n) quantile(params$f_predict[,n], probs=probs))

  plot(1, type="n", main=title,
       xlab="x", ylab="y", xlim=c(-10, 10), ylim=c(-10, 10))
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
  lines(true$x, true$f, lwd=2, xlab="x", ylab="y", col="black")
  points(data$x, data$y, col="white", pch=16, cex=1.2)
  points(data$x, data$y, col="black", pch=16, cex=0.8)
}

# Plot Gaussian process predictive quantiles
plot_gp_pred_quantiles <- function(fit, data, true, title) {
  params <- extract(fit)

  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:length(data$x_predict),
                 function(n) quantile(params$y_predict[,n], probs=probs))

  plot(1, type="n", main=title,
       xlab="x", ylab="y", xlim=c(-10, 10), ylim=c(-10, 10))
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
  lines(true$x, true$f, lwd=2, xlab="x", ylab="y", col="black")
  points(data$x, data$y, col="white", pch=16, cex=1.2)
  points(data$x, data$y, col="black", pch=16, cex=0.8)
}
