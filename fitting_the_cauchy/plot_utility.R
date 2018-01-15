c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

breaks <- function(mag, N) {
  return(mag*(2*(0:N)/N - 1))
}

# Plot posterior quantiles
plot_estimated_quantiles <- function(fit, title) {
  params <- extract(fit)

  probs = c(0.05, 0.5, 0.95)
  cred <- sapply(1:50, function(n) quantile(params$x[,n], probs=probs))

  p1 <- hist(cred[1,], breaks=breaks(10, 100), plot=FALSE)
  p2 <- hist(cred[2,], breaks=breaks(10, 100), plot=FALSE)
  p3 <- hist(cred[3,], breaks=breaks(10, 100), plot=FALSE)

  plot(p1, col=c_dark, border=c_dark_highlight, main=title, xlab="x", yaxt='n', ylim=c(0, 32), ylab="")
  plot(p2, col=c_dark, border=c_dark_highlight, add=T)
  plot(p3, col=c_dark, border=c_dark_highlight, add=T)

  abline(v=qcauchy(0.05, location = 0, scale = 1), col=c_light, lty=1, lwd=2)
  abline(v=qcauchy(0.50, location = 0, scale = 1), col=c_light, lty=1, lwd=2)
  abline(v=qcauchy(0.95, location = 0, scale = 1), col=c_light, lty=1, lwd=2)

  rect(-8, 31, 8, 35, border=NA, col="white")

  text(x=qcauchy(0.05, location = 0, scale = 1), y=32, labels="5%")
  text(x=qcauchy(0.50, location = 0, scale = 1), y=32, labels="50%")
  text(x=qcauchy(0.95, location = 0, scale = 1), y=32, labels="95%")
}
