data {
  int<lower=1> N;      // Number of observations
  real y[N];           // Observations
  real<lower=0> sigma; // Measurement variability
}

parameters {
  real a;
}

model {
  // Prior model
  a ~ normal(0, 100);
  
  // Observational model
  y ~ normal(a, sigma);
}
