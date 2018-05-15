data {
  int N;
}

generated quantities {
  real<lower=0> lambda = fabs(normal_rng(0, 6.44787));
  int y[N];
    for (n in 1:N) y[n] = poisson_rng(lambda);
}
