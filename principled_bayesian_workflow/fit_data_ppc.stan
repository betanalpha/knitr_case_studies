data {
  int N;    // Number of observations
  int y[N]; // Count at each observation
}

parameters {
  real<lower=0> lambda; // Poisson intensity
}

model {
  lambda ~ normal(0, 6.44787); // Prior model
  y ~ poisson(lambda);         // Observational model
}

// Simulate a full observation from the current value of the parameters
generated quantities {
  int y_ppc[N];

  for (n in 1:N)
    y_ppc[n] = poisson_rng(lambda);
}
