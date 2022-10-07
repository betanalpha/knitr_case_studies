data {
  int<lower=0> K;                       // Number of contexts
  int<lower=0> N;                       // Number of observations
  vector[N] y;                          // Observations
  int<lower=1, upper=K> context_idx[N]; // Context assignments
  real<lower=0> sigma;                  // Measurement variability
}

parameters {
  // Laplace parameters
  vector[K] theta;
}

model {
  // Laplace prior model
  theta ~ double_exponential(0, 0.1);

  // Observational model
  y ~ normal(theta[context_idx], sigma);
}
