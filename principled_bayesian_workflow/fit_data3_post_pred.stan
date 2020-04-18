data {
  int N;    // Number of observations
  int y[N]; // Count at each observation
}

parameters {
  real<lower=0> lambda;         // Poisson intensity
  real<lower=0, upper=1> theta; // Excess zero probability
}

model {
  // Prior model
  lambda ~ inv_gamma(3.48681, 9.21604);
  theta ~ beta(2.8663, 2.8663);

  // Observational model that mixes a Poisson with excess zeros
  for (n in 1:N) {
    real lpdf = poisson_lpmf(y[n] | lambda);
    if (y[n] == 0)
      target += log_mix(theta, 0, lpdf);
    else
      target += log(1 - theta) + lpdf;
  }
}

// Simulate a full observation from the current value of the parameters
generated quantities {
  int y_post_pred[N] = rep_array(0, N);

  for (n in 1:N)
    if (!bernoulli_rng(theta))
      y_post_pred[n] = poisson_rng(lambda);
}
