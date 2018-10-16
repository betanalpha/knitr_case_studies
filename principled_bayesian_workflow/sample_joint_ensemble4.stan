data {
  int N;
}

transformed data {
  int U = 14;
}

generated quantities {
  // Simulate model configuration from prior model
  real<lower=0, upper=1> theta = beta_rng(2.8663, 2.8663);
  real<lower=0> lambda = inv_gamma_rng(3.48681, 9.21604);

  // Simulate data from observational model
  int y[N] = rep_array(0, N);
  for (n in 1:N) {
    if (!bernoulli_rng(theta)) {
      real sum_p = 0;
      real u = uniform_rng(0, 1);

      for (b in 0:U) {
        sum_p = sum_p + exp(poisson_lpmf(b | lambda) - poisson_lcdf(U | lambda));
        if (sum_p >= u) {
          y_ppc[n] = b;
          break;
        }
      }
    }
  }
}
