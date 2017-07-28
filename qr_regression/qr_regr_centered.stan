data {
  int<lower=1> N;
  int<lower=1> M;
  matrix[M, N] X;
  vector[N] y;
}

transformed data {
  matrix[M, N] X_centered;
  matrix[N, M] Q;
  matrix[M, M] R;
  matrix[M, M] R_inv;

  for (m in 1:M)
    X_centered[m] = X[m] - mean(X[m]);

  // Compute, thin, and then scale QR decomposition
   Q = qr_Q(X_centered')[, 1:M] * N;
   R = qr_R(X_centered')[1:M, ] / N;
   R_inv = inverse(R);
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
  alpha ~ normal(0, 100);
  sigma ~ cauchy(0, 10);

  y ~ normal(Q * beta_tilde + alpha, sigma);
}
