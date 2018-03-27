g <- function(x, a) 10 * exp(-0.1 * (x - 2.5)**2 ) + 8 * exp(-0.2 * (x + 4)**2) - a;

middle_min_g <- 4.211305
middle_max_g <- 8.152303
max_g <- 10.00171


#heights <- c(5)
#heights <- seq(2, 10, 2)
heights <- seq(0.2, 10, 0.2)

height <- c()
lower <- c()
upper <- c()

for (y in heights) {
  if (y < middle_min_g) {
    height <- c(height, y)
    lower <- c(lower, uniroot(function(x) g(x, y), c(-10, -4), maxiter = 10000)$root)
    upper <- c(upper, uniroot(function(x) g(x, y), c(2.5, 10), maxiter = 10000)$root)
  } else {
    if (y < middle_max_g) {
      height <- c(height, y)
      lower <- c(lower, uniroot(function(x) g(x, y), c(-10, -4), maxiter = 10000)$root)
      upper <- c(upper, uniroot(function(x) g(x, y), c(-4, -1.2), maxiter = 10000)$root)
    }
    
    height <- c(height, y)
    lower <- c(lower, uniroot(function(x) g(x, y), c(-1.2, 2.5), maxiter = 10000)$root)
    upper <- c(upper, uniroot(function(x) g(x, y), c(2.5, 10), maxiter = 10000)$root)
  }
}

toString(rev(lower))
toString(rev(upper))
toString(rev(height))
