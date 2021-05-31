# Plot expected daily sales
plot_excess_variation <- function(samples, data, truth, title) {

  idx <- rep(1:data$N_obs_days, each=2)
  x <- sapply(1:length(idx), function(k) if(k %% 2 == 0) idx[k] + 0.5 else idx[k] - 0.5)

  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:data$N_obs_days, function(k) quantile(samples$delta_days[,k], probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))

  plot(1, type="n", main=title,
       xlim=c(0.5, data$N_obs_days + 0.5), xlab="Day",
       ylim=c(-0.5, 0.5), ylab="Inferred Excess Sales Variation")

  polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  for (k in 1:data$N_obs_days) {
    lines(x[(2 * k - 1):(2 * k)], pad_cred[5,(2 * k - 1):(2 * k)], col=c_dark, lwd=2)
  }
  
  points(1:data$N_obs_days, truth$delta_days, col=c_mid_teal, pch=16, cex=0.5)
}

# Plot expected daily sales residual
plot_excess_variation_residual <- function(samples, data, truth, title) {
  
  idx <- rep(1:data$N_obs_days, each=2)
  x <- sapply(1:length(idx), function(k) if(k %% 2 == 0) idx[k] + 0.5 else idx[k] - 0.5)
  
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:data$N_obs_days, function(k) quantile(samples$delta_days[,k] - truth$delta_days[k], probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))
  
  plot(1, type="n", main=title,
       xlim=c(0.5, data$N_obs_days + 0.5), xlab="Day",
       ylim=c(-0.5, 0.5), ylab="Inferred Excess Sales Variation Minus Truth")
  
  abline(h=0, col="gray80", lty=1, lw=2)
  
  polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  for (k in 1:data$N_obs_days) {
    lines(x[(2 * k - 1):(2 * k)], pad_cred[5,(2 * k - 1):(2 * k)], col=c_dark, lwd=2)
  }
  
}

# Plot expected daily sales
plot_expected_sales <- function(samples, data, title) {

  idx <- rep(1:data$N_obs_days, each=2)
  x <- sapply(1:length(idx), function(k) if(k %% 2 == 0) idx[k] + 0.5 else idx[k] - 0.5)

  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:data$N_obs_days, function(k) quantile(exp(samples$log_mu[,k]), probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))

  plot(1, type="n", main=title,
       xlim=c(0.5, data$N_obs_days + 0.5), xlab="Day",
       ylim=c(0, 400), ylab="Expected Sales")

  polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  for (k in 1:data$N_obs_days) {
    lines(x[(2 * k - 1):(2 * k)], pad_cred[5,(2 * k - 1):(2 * k)], col=c_dark, lwd=2)
  }

  points(data$obs_day_idx, data$y_obs, col="black", pch=16, cex=0.25)
}

# Plot residual expected daily sales
plot_expected_sales_residual <- function(samples, data, title) {

  idx <- rep(1:data$N_obs_days, each=2)
  x <- sapply(1:length(idx), function(k) if(k %% 2 == 0) idx[k] + 0.5 else idx[k] - 0.5)
  
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:data$N_obs_days, function(k) quantile(exp(samples$log_mu[,k]), probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))
  
  plot(1, type="n", main=title,
       xlim=c(0.5, data$N_obs_days + 0.5), xlab="Day",
       ylim=c(-100, 100), ylab="Expected Sales Residuals")
  
  polygon(c(x, rev(x)), c(pad_cred[1,] - pad_cred[5,], rev(pad_cred[9,] - pad_cred[5,])),
          col = c_light, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[2,] - pad_cred[5,], rev(pad_cred[8,] - pad_cred[5,])),
          col = c_light_highlight, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[3,] - pad_cred[5,], rev(pad_cred[7,] - pad_cred[5,])),
          col = c_mid, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[4,] - pad_cred[5,], rev(pad_cred[6,] - pad_cred[5,])),
          col = c_mid_highlight, border = NA)
  for (k in 1:data$N_obs_days) {
    lines(x[(2 * k - 1):(2 * k)], c(0, 0), col=c_dark, lwd=2)
  }
  
  data_residual <- sapply(1:data$N_obs, function(n) data$y_obs[n] - cred[5, data$obs_day_idx[n]])
  points(data$obs_day_idx, data_residual, col="black", pch=16, cex=0.25)
}
  
# Plot retrodictions and predictions
plot_dictions <- function(samples, data, title) {

  idx <- rep(1:data$N_obs_days, each=2)
  x <- sapply(1:length(idx), function(k) if(k %% 2 == 0) idx[k] + 0.5 else idx[k] - 0.5)
  
  plot(1, type="n", main=title,
       xlim=c(0.5, 500.5), xlab="Day",
       ylim=c(0, 400), ylab="Retrodicted/Predicted Sales")
  
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:data$N_obs_days, function(k) quantile(samples$y_post_retro[,k], probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))
  
  polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  for (k in 1:data$N_obs_days) {
    lines(x[(2 * k - 1):(2 * k)], pad_cred[5,(2 * k - 1):(2 * k)], col=c_dark, lwd=2)
  }
  
  points(data$obs_day_idx, data$y_obs, col="black", pch=16, cex=0.25)
  
  idx <- rep(data$pred_day_idx, each=2)
  x <- sapply(1:length(idx), function(k) if(k %% 2 == 0) idx[k] + 0.5 else idx[k] - 0.5)
  
  cred <- sapply(1:data$N_pred_days, function(k) quantile(samples$y_post_pre[,k], probs=probs))
  pad_cred <- do.call(cbind, lapply(idx - idx[1] + 1, function(n) cred[1:9,n]))
  
  polygon(c(x, rev(x)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  for (k in 1:data$N_pred_days) {
    lines(x[(2 * k - 1):(2 * k)], pad_cred[5,(2 * k - 1):(2 * k)], col=c_dark, lwd=2)
  }
  
  points(data$pred_day_idx, data$y_pred, col=c_mid_teal, pch=16, cex=0.25)
  
  text(10, 50, cex=1.25, label="Posterior\nRetrodictions", 
       pos=4, col="black")
  text(data$pred_day_idx[1] + 0, 50, cex=1.25, label="Posterior\nPredictions", 
       pos=4, col=c_mid_teal)
}
    
# Plot retrodiction and prediction residuals
plot_dictions_residual <- function(samples, data, title) {
  
  idx <- rep(1:data$N_obs_days, each=2)
  x <- sapply(1:length(idx), function(k) if(k %% 2 == 0) idx[k] + 0.5 else idx[k] - 0.5)
  
  plot(1, type="n", main=title,
       xlim=c(0.5, 500.5), xlab="Day",
       ylim=c(-100, 100), ylab="Retrodicted/Predicted Sales Residuals")
  
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:data$N_obs_days, function(k) quantile(samples$y_post_retro[,k], probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))
  
  polygon(c(x, rev(x)), c(pad_cred[1,] - pad_cred[5,], rev(pad_cred[9,] - pad_cred[5,])),
          col = c_light, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[2,] - pad_cred[5,], rev(pad_cred[8,] - pad_cred[5,])),
          col = c_light_highlight, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[3,] - pad_cred[5,], rev(pad_cred[7,] - pad_cred[5,])),
          col = c_mid, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[4,] - pad_cred[5,], rev(pad_cred[6,] - pad_cred[5,])),
          col = c_mid_highlight, border = NA)
  for (k in 1:data$N_obs_days) {
    lines(x[(2 * k - 1):(2 * k)], c(0, 0), col=c_dark, lwd=2)
  }
  
  data_residual <- sapply(1:data$N_obs, function(n) data$y_obs[n] - cred[5, data$obs_day_idx[n]])
  points(data$obs_day_idx, data_residual, col="black", pch=16, cex=0.25)
  
  idx <- rep(data$pred_day_idx, each=2)
  x <- sapply(1:length(idx), function(k) if(k %% 2 == 0) idx[k] + 0.5 else idx[k] - 0.5)
  
  cred <- sapply(1:data$N_pred_days, function(k) quantile(samples$y_post_pre[,k], probs=probs))
  pad_cred <- do.call(cbind, lapply(idx - idx[1] + 1, function(n) cred[1:9,n]))
  
  polygon(c(x, rev(x)), c(pad_cred[1,] - pad_cred[5,], rev(pad_cred[9,] - pad_cred[5,])),
          col = c_light, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[2,] - pad_cred[5,], rev(pad_cred[8,] - pad_cred[5,])),
          col = c_light_highlight, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[3,] - pad_cred[5,], rev(pad_cred[7,] - pad_cred[5,])),
          col = c_mid, border = NA)
  polygon(c(x, rev(x)), c(pad_cred[4,] - pad_cred[5,], rev(pad_cred[6,] - pad_cred[5,])),
          col = c_mid_highlight, border = NA)
  for (k in 1:data$N_pred_days) {
    lines(x[(2 * k - 1):(2 * k)], c(0, 0), col=c_dark, lwd=2)
  }
  
  points(data$pred_day_idx, data$y_pred - cred[5,], col=c_mid_teal, pch=16, cex=0.25)
  
  text(10, -75, cex=1.25, label="Posterior\nRetrodictions", 
       pos=4, col="black")
  text(data$pred_day_idx[1] + 0, -75, cex=1.25, label="Posterior\nPredictions", 
       pos=4, col=c_mid_teal)
}

# Empirical continuous rank probability score
empirical_crps <- function(xs, x_tilde) {
  
  xs_sorted <- sort(xs)
  N <- length(xs_sorted)
  n_tilde <- sum(xs_sorted < x_tilde) + 1
  
  if (n_tilde == 1) {
    # All forecasts are above x
    xs_expanded <- c(x_tilde, xs_sorted)
    deltas <- sapply(1:N, function(n) xs_expanded[n + 1] - xs_expanded[n])
    ps <- (0:(N - 1)) / N
    rps <- sum( (1 - ps[n_tilde:N])**2 * deltas[n_tilde:N])
  } else if (n_tilde == N + 1) {
    # All forecasts are below x
    xs_expanded <- c(xs_sorted, x_tilde)
    deltas <- sapply(1:N, function(n) xs_expanded[n + 1] - xs_expanded[n])
    ps <- (1:N) / N
    rps <- sum(ps[1:(n_tilde - 1)]**2 * deltas[1:(n_tilde - 1)]) 
  } else {
    # Forecasts surround x
    xs_expanded <- c(xs_sorted[1:(n_tilde - 1)], x_tilde, xs_sorted[n_tilde:N])
    deltas <- sapply(1:N, function(n) xs_expanded[n + 1] - xs_expanded[n])
    ps <- c(1:(n_tilde - 1), n_tilde:N - 1) / N
    rps <- sum(ps[1:(n_tilde - 1)]**2 * deltas[1:(n_tilde - 1)]) +
           sum( (1 - ps[n_tilde:N])**2 * deltas[n_tilde:N])
  } 
  
  return(rps)
}

# Summed continuous rank probability score
sum_empirical_crps <- function(x_forecast, x_truth) {
  sum <- 0
  for (k in 1:length(x_truth)) {
    xs <- x_forecast[,k]
    x_tilde <- x_truth[k]
    sum <- sum + empirical_crps(xs, x_tilde)
  }
  return(sum)
}
 
