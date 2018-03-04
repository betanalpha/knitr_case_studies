c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

# Plot posterior quantiles
plot_post_quantiles <- function(fit, data, title) {
  large_slope_idx = which(abs(input_data$beta_true) > 3)

  params = extract(fit)

  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:input_data$M, function(m) quantile(params$beta[,m], probs=probs))

  plot(1, type="n", main=title,
       xlim=c(1, input_data$M), xlab="Slope Index",
       ylim=c(min(cred[1,]), max(cred[9,])),
       ylab="Marginal Slope Posteriors")
  sapply(large_slope_idx, function(idx) abline(v=idx, col="gray80", lwd=2, lty=3))

  polygon(c(1:input_data$M, input_data$M:1), c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(1:input_data$M, input_data$M:1), c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(1:input_data$M, input_data$M:1), c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(1:input_data$M, input_data$M:1), c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(1:input_data$M, cred[5,], col=c_dark, lwd=2)

  lines(1:input_data$M, input_data$beta_true, lwd=1.5, col="white")
  lines(1:input_data$M, input_data$beta_true, lwd=1.25, col="black")
}

# Plot residual quantiles
plot_residual_quantiles <- function(fit, data, title) {
  large_slope_idx = which(abs(input_data$beta_true) > 3)

  params = extract(fit)

  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:input_data$M, function(m) quantile(params$beta[,m] - input_data$beta_true[m], probs=probs))

  plot(1, type="n", main=title,
       xlim=c(1, input_data$M), xlab="Slope Index",
       ylim=c(min(cred[1,]), max(cred[9,])),
       ylab="Marginal Slope Posteriors Relative to Truth")
  sapply(large_slope_idx, function(idx) abline(v=idx, col="gray80", lwd=2, lty=3))

  polygon(c(1:input_data$M, input_data$M:1), c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(1:input_data$M, input_data$M:1), c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(1:input_data$M, input_data$M:1), c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(1:input_data$M, input_data$M:1), c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(1:input_data$M, cred[5,], col=c_dark, lwd=2)

  abline(h=0, col="white", lwd=1.5)
  abline(h=0, col="black", lwd=1.25)
}

# Plot auxiliary marginal posteriors
plot_aux_posteriors <- function(fit, title) {
  params = extract(fit)

  par(mfrow=c(1, 2))

  hist(params$alpha, main=title, col=c_dark, border=c_dark_highlight,
       xlab="alpha", yaxt='n', ylab="")
  abline(v=3,col="white", lty=1, lwd=2.5)
  abline(v=3,col="black", lty=1, lwd=2)

  hist(params$sigma, main=title, col=c_dark, border=c_dark_highlight,
       xlab="sigma", yaxt='n', ylab="")
  abline(v=1,col="white", lty=1, lwd=2.5)
  abline(v=1,col="black", lty=1, lwd=2)
}
