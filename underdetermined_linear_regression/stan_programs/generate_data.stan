transformed data {
  int<lower=0> M = 200;
  int<lower=0> N = 100;
  real alpha = 3;
  real sigma = 1;
}

generated quantities {
  matrix[N, M] X;
  real y[N];

  vector[M] beta;
  for (m in 1:M) beta[m] = normal_rng(0, 10);

  for (n in 1:N) {
    for (m in 1:M) X[n, m] = normal_rng(0, 1);
    y[n] = normal_rng(X[n, 1:M] * beta + alpha, sigma);
  }
}
