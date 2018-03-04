data {
  int<lower=1> N; // Number of data
  int<lower=1> M; // Number of covariates
  matrix[M, N] X;
  real y[N];
}

parameters {
  vector[M] beta;
  real alpha;
  real<lower=0> sigma;
}

model {
  // Weakly regularizing priors on the slopes
  beta ~ normal(0, 10);
  alpha ~ normal(0, 2);
  sigma ~ normal(0, 2);

  y ~ normal(X' * beta + alpha, sigma);
}
