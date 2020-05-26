data {
  int<lower=1> N;  // Number of observations
  real y[N];       // Observations
  
  int<lower=1> N_aux;  // Number of auxiliary observations
  real y_aux[N_aux];   // Auxiliary observations
}

parameters {
  real a;
  real b;
  real<lower=0> sigma;
  real<lower=0> sigma_aux;
}

model {
  // Prior model
  a ~ normal(0, 1);
  b ~ normal(0, 1);
  sigma ~ normal(0, 1);
  sigma_aux ~ normal(0, 1);
  
  // Auxiliary observational model
  y_aux ~ normal(a, sigma_aux);
  
  // Observational model
  y ~ normal(a + b, sigma);
}
