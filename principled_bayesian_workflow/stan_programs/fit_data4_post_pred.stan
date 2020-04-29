data {
  int N;    // Number of observations
  int y[N]; // Count at each observation
}

transformed data {
  int U = 14;  // Observational cutoff
}

parameters {
  real<lower=0> lambda;         // Poisson intensity
  real<lower=0, upper=1> theta; // Excess zero probability
}

model {
  // Prior model
  lambda ~ inv_gamma(3.5, 9);
  theta ~ beta(3, 3);

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
  int y_post_pred[N] = rep_array(0, N);

  for (n in 1:N) {
    if (!bernoulli_rng(theta)) {
      real sum_p = 0;
      real u = uniform_rng(0, 1);

      for (b in 0:U) {
        sum_p = sum_p + exp(poisson_lpmf(b | lambda) - poisson_lcdf(U | lambda));
        if (sum_p >= u) {
          y_post_pred[n] = b;
          break;
        }
      }
    }
  }
}
