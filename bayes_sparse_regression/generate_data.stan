transformed data {
  int<lower=0> M = 200;
  int<lower=0> N = 100;
  real alpha = 3;
  real sigma = 1;
  real<lower=0, upper=1> sig_prob = 0.05;
}

generated quantities {
  matrix[M, N] X;
  real y[N];

  vector[M] beta;
  for (m in 1:M) {
    if (bernoulli_rng(sig_prob))
      if (bernoulli_rng(0.5))
        beta[m] = normal_rng(10, 1);
      else
        beta[m] = normal_rng(-10, 1);
    else
      beta[m] = normal_rng(0, 0.25);
  }

  for (n in 1:N) {
    for (m in 1:M)
      X[m, n] = normal_rng(0, 1);

    y[n] = normal_rng(X[1:M,n]' * beta + alpha, sigma);
  }
}
