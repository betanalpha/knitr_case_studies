data {
  int<lower=1> N; // Number of observations
  int<lower=1> M; // Number of covariates
  matrix[N, M] X; // Design matrix
  vector[N] y;    // Observations
}

parameters {
  vector[M] beta;
  real alpha;
  real<lower=0> sigma;
}

model {
  // Weakly informative containment priors
  beta ~ normal(0, 10);
  alpha ~ normal(0, 2);
  sigma ~ normal(0, 2);

  // Observational model
  y ~ normal(X * beta + alpha, sigma);
}
