data {
  int<lower=0> K;
}

parameters {
  real<lower=0> tau;
  real phi[K];
}

model {
 tau ~ normal(0, 5);
 phi ~ normal(0, tau);
}
