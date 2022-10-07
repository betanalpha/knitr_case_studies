data {
  int<lower=0> K;                       // Number of contexts
  int<lower=0> N;                       // Number of observations
  vector[N] y;                          // Observations
  int<lower=1, upper=K> context_idx[N]; // Context assignments
  real<lower=0> sigma;                  // Measurement variability
}

parameters {
  // Cauchy parameters
  vector[K] theta;
}

model {
  // Cauchy prior model
  theta ~ cauchy(0, 0.1);

  // Observational model
  y ~ normal(theta[context_idx], sigma);
}
