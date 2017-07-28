data {
  int<lower=1> N;
  int<lower=1> M;
  matrix[M, N] X;
  vector[N] y;
}

transformed data {
  // Compute, thin, and then scale QR decomposition
  matrix[N, M] Q = qr_Q(X')[, 1:M] * N;
  matrix[M, M] R = qr_R(X')[1:M, ] / N;
  matrix[M, M] R_inv = inverse(R);
}

parameters {
  vector[M] beta_tilde;
  real alpha;
  real<lower=0> sigma;
}

transformed parameters {
  vector[M] beta = R_inv * beta_tilde;
}

model {
  beta ~ normal(0, 10);
  alpha ~ normal(0, 10);
  sigma ~ cauchy(0, 10);

  y ~ normal(Q * beta_tilde + alpha, sigma);
}
