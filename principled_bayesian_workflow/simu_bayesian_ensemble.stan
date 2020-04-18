data {
  int N;
}

generated quantities {
  // Simulate model configuration from prior model
  real<lower=0> lambda = fabs(normal_rng(0, 5.82337));
  
  // Simulate data from observational model
  int y[N];
  for (n in 1:N) y[n] = poisson_rng(lambda);
}
