data {
  int N;
}

generated quantities {
  // Simulate model configuration from prior model
  real<lower=0> lambda = inv_gamma_rng(3.5, 9);
  real<lower=0, upper=1> theta = beta_rng(3, 3);
  
  // Simulate data from observational model
  int y[N] = rep_array(0, N);
  for (n in 1:N)
    if (!bernoulli_rng(theta))
      y[n] = poisson_rng(lambda);
}
