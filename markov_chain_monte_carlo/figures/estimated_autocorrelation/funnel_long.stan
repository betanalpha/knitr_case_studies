parameters {
  real x;
  real<lower=0> tau;
}

model {
  x ~ normal(0, tau);
  tau ~ normal(0, 3);
}
