data {
  int<lower=1> N;
  int<lower=1> M;
  matrix[M, N] X;
  vector[N] y;
}

parameters {
  vector[M] beta;
  real alpha;
  real<lower=0> sigma;
}

model {
  beta ~ normal(0, 1);
  alpha ~ normal(0, 1);
  sigma ~ cauchy(0, 1);

  y ~ normal(X' * beta + alpha, sigma);
}
