data {
  int<lower=1> N;
  int<lower=1> M;
  matrix[N, M] Q;
  matrix[M, M] R;
  vector[N] y;
}

parameters {
  vector[M] beta_tilde;
  real alpha;
  real<lower=0> sigma;
}

transformed parameters {
  // Lots of transposing because Stan doesn't
  // have a mdivide_right_tri_upper
  vector[M] beta = mdivide_right_tri_low(beta_tilde', R')';
}

model {
  beta_tilde ~ normal(0, 10);
  alpha ~ normal(0, 10);
  sigma ~ cauchy(0, 10);

  y ~ normal(Q * beta_tilde + alpha, sigma);
}
