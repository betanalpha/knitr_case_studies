data {
  int N;
}

transformed data {
  real lambda = 9;
  real theta = 0.3;
  int U = 14;
}

generated quantities {
  int y[N] = rep_array(0, N);
  for (n in 1:N) {
    if (!bernoulli_rng(theta)) {
      y[n] = U + 1;
      while (y[n] > U)
        y[n] = poisson_rng(lambda);
    }
  }
}
