data {
  int<lower=1> N;      // Number of observations
  real y[N];           // Observations
  real<lower=0> sigma; // Measurement variability
  
  int<lower=1> N_aux;      // Number of auxiliary observations
  real y_aux[N_aux];       // Auxiliary observations
  real<lower=0> sigma_aux; // Auxiliary measurement variability
}

parameters {
  real a;
  real b;
}

model {
  // Prior model
  a ~ normal(0, 100);
  b ~ normal(0, 100);
  
  // Auxiliary observational model
  y_aux ~ normal(a, sigma_aux);
  
  // Observational model
  y ~ normal(a + b, sigma);
}
