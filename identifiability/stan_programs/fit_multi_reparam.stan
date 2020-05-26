data {
  int<lower=1> N;      // Number of observations
  real y[N];           // Observations
  real<lower=0> sigma; // Measurement variability
}

parameters {
  real a_times_b;
}

model {
  // Prior model
  a_times_b ~ normal(0, 100);
  
  // Observational model
  y ~ normal(a_times_b, sigma);
}
