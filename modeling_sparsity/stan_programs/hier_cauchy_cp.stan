data {
  int<lower=0> K;                       // Number of contexts
  int<lower=0> N;                       // Number of observations
  vector[N] y;                          // Observations
  int<lower=1, upper=K> context_idx[N]; // Context assignments
  real<lower=0> sigma;                  // Measurement variability
}

parameters {
  // Horseshoe parameters
  vector[K] theta;
  real<lower=0> tau;
}

model {
  // Horseshoe prior model
  theta ~ cauchy(0, tau);
  tau ~ normal(0, 10);

  // Observational model
  y ~ normal(theta[context_idx], sigma);
}
