data {
  int N;
}

generated quantities {
  real<lower=0, upper=1> theta = beta_rng(1, 1);
  real<lower=0> lambda = fabs(normal_rng(0, 6.44787));

  int y[N] = rep_array(0, N);
  for (n in 1:N)
    if (!bernoulli_rng(theta))
      y[n] = poisson_rng(lambda);
}
