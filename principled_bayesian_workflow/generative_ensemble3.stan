data {
  int N;
}

generated quantities {
  real<lower=0, upper=1> theta = beta_rng(2.8663, 2.8663);
  real<lower=0> lambda = inv_gamma_rng(3.48681, 9.21604);

  int y[N] = rep_array(0, N);
  for (n in 1:N)
    if (!bernoulli_rng(theta))
      y[n] = poisson_rng(lambda);
}
