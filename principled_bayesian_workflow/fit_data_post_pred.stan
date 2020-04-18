data {
  int N;    // Number of observations
  int y[N]; // Count at each observation
}

parameters {
  real<lower=0> lambda; // Poisson intensity
}

model {
  lambda ~ normal(0, 5.82337); // Prior model
  y ~ poisson(lambda);         // Observational model
}

// Simulate a full observation from the current value of the parameters
generated quantities {
  int y_post_pred[N];
  for (n in 1:N)
    y_post_pred[n] = poisson_rng(lambda);
}
