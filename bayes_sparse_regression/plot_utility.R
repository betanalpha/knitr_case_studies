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

  idx <- rep(1:input_data$M, each=2)
  x <- sapply(1:length(idx), function(m) if(m %% 2 == 0) idx[m] + 0.5 else idx[m] - 0.5)
  pad_beta_true <- do.call(cbind, lapply(idx, function(n) beta_true[n]))

  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:input_data$M, function(m) quantile(params$beta[,m], probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))

min(c(beta_true, cred[1,]))

  plot(1, type="n", main=title,
       xlim=c(0.5, input_data$M + 0.5), xlab="Slope Index",
       ylim=c(min(c(beta_true, cred[1,])), max(c(beta_true, cred[9,]))),
       ylab="Slope Posterior")
  sapply(large_slope_idx, function(idx) abline(v=idx, col="gray80", lwd=2, lty=3))

  polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(x, pad_cred[5,], col=c_dark, lwd=2)

  lines(x, pad_beta_true, lwd=1.5, col="white")
  lines(x, pad_beta_true, lwd=1.25, col="black")
}

# Plot residual quantiles
plot_residual_quantiles <- function(fit, data, title) {
  large_slope_idx = which(abs(input_data$beta_true) > 3)

  params = extract(fit)

  idx <- rep(1:input_data$M, each=2)
  x <- sapply(1:length(idx), function(m) if(m %% 2 == 0) idx[m] + 0.5 else idx[m] - 0.5)
  pad_beta_true <- do.call(cbind, lapply(idx, function(n) beta_true[n]))

  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:input_data$M, function(m) quantile(params$beta[,m] - input_data$beta_true[m], probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))

  plot(1, type="n", main=title,
       xlim=c(0.5, input_data$M + 0.5), xlab="Slope Index",
       ylim=c(min(c(beta_true, cred[1,])), max(c(beta_true, cred[9,]))),
       ylab="Slope Posterior Minus True Slope")
  sapply(large_slope_idx, function(idx) abline(v=idx, col="gray80", lwd=2, lty=3))

  polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(x, pad_cred[5,], col=c_dark, lwd=2)

  abline(h=0, col="white", lwd=1.5)
  abline(h=0, col="black", lwd=1.25)
}

# Plot summary residual quantiles
plot_summary_residual_quantiles <- function(fit, data, title) {
  large_slope_idx = which(abs(input_data$beta_true) > 3)

  params = extract(fit)

  idx <- rep(51:100, each=2)
  x <- sapply(1:length(idx), function(m) if(m %% 2 == 0) idx[m] + 0.5 else idx[m] - 0.5)
  pad_beta_true <- do.call(cbind, lapply(idx, function(n) beta_true[n]))

  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:input_data$M, function(m) quantile(params$beta[,m] - input_data$beta_true[m], probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))

  plot(1, type="n", main=title,
       xlim=c(50.5, 100.5), xlab="Slope Index",
       ylim=c(-3, 3),
       ylab="Slope Posterior Minus True Slope")
  sapply(large_slope_idx, function(idx) abline(v=idx, col="gray80", lwd=2, lty=3))

  polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(x, pad_cred[5,], col=c_dark, lwd=2)

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
