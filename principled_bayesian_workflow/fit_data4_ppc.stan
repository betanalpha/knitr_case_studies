data {
  int N;    // Number of observations
  int y[N]; // Count at each observation
}

transformed data {
  int U = 14;  // Observational cutoff
}

parameters {
  real<lower=0, upper=1> theta; // Excess zero probability
  real<lower=0> lambda;         // Poisson intensity
}

model {
  // Prior model
  theta ~ beta(2.8663, 2.8663);
  lambda ~ inv_gamma(3.48681, 9.21604);

  // Observational model that mixes a truncated Poisson with excess zeros
  for (n in 1:N) {
    real lpdf = poisson_lpmf(y[n] | lambda) - poisson_lcdf(U | lambda);
    if (y[n] == 0)
      target += log_mix(theta, 0, lpdf);
    else
      target += log(1 - theta) + lpdf;
  }
}

// Simulate a full observation from the current value of the parameters
generated quantities {
  int y_ppc[N] = rep_array(0, N);

  for (n in 1:N) {
    if (!bernoulli_rng(theta)) {
      y_ppc[n] = U + 1;
      while (y_ppc[n] > U)
        y_ppc[n] = poisson_rng(lambda);
    }
  }
}
