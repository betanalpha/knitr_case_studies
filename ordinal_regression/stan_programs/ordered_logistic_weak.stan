data {
  int<lower=1> N;             // Number of observations
  int<lower=1> K;             // Number of ordinal categories
  int<lower=1, upper=K> y[N]; // Observed ordinals
}

parameters {
  real gamma;       // Latent effect
  ordered[K - 1] c; // (Internal) cut points
}

model {
  // Prior model
  gamma ~ normal(0, 1);
  
  // Observational model
  y ~ ordered_logistic(rep_vector(gamma, N), c);
}

generated quantities {
  int<lower=1, upper=K> y_ppc[N];
  for (n in 1:N)
    y_ppc[n] = ordered_logistic_rng(gamma, c);
}
