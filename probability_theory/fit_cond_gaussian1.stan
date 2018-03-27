transformed data {
  real x1 = 1;
  real mu1 = 2;
  real sigma1 = 0.75;

  real mu2 = -1;
  real sigma2 = 1.2;

  real rho = 0.7;

  real cond_mu = mu2 + (sigma2 / sigma1) * rho * (x1 - mu1);
  real cond_sd = sigma2 * sqrt(1 - square(rho));
}

parameters {
  real x2;
}

model {
  x2 ~ normal(cond_mu, cond_sd);
}
