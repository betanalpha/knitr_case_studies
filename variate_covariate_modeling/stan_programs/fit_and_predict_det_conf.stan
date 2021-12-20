functions {
  // Location-dispersion paramaeterization of gamma family
  real gamma_ld_lpdf(vector y, vector mu, real phi) {
    real beta = inv(phi);
    vector[num_elements(mu)] alpha = beta * mu;
    return gamma_lpdf(y | alpha, beta);
  }
  
  // Location-dispersion paramaeterization of gamma family
  real gamma_ld_rng(real mu, real phi) {
    real beta = inv(phi);
    real alpha = beta * mu;
    return gamma_rng(alpha, beta);
  }
}

data {
  int<lower=1> N1; // Number of complete observations
  vector[N1] x1;   // Complete observations
  vector[N1] y1;   // Complete observations
  
  int<lower=1> N2; // Number of incomplete observations
  vector[N2] x2;   // Incomplete observations
  
  int<lower=1> N_grid; // Number of grid points for quantifying functional behavior
  vector[N_grid] x_grid; // Grid points for quantifying functional behavior
}

parameters { 
  real<lower=0> lambda; // Detector calibration, output units / input units
  real<lower=0> psi;    // Detector resolution

  real<lower=0> gamma;  // Detector resolution correction factor (inverse input units)
  
  // Location and dispersion of completely observed covariate distribution 
  // Location has input units, dispersion is unitless
  real<lower=0> mu1;
  real<lower=0> phi1;
  
  // Location and dispersion of incompletely observed covariate distribution
  // Location has input units, dispersion is unitless
  real<lower=0> mu2; 
  real<lower=0> phi2;
}

model {
  lambda ~ normal(0, 0.2);
  psi ~ normal(0, 0.2);
  
  gamma ~ normal(0, 10);
  
  mu1 ~ normal(0, 50);
  phi1 ~ normal(0, 5);
  
  mu2 ~ normal(0, 50);
  phi2 ~ normal(0, 5);

  x1 ~ gamma_ld(rep_vector(mu1, N1), phi1);
  y1 ~ gamma_ld(lambda * x1, psi + 0.001 * gamma * mu1);
  
  x2 ~ gamma_ld(rep_vector(mu2, N2), phi2);
}

generated quantities {
  real x1_pred[N1]; // Complete observation covariate retrodiction
  real y1_pred[N1]; // Complete observation conditional variate retrodiction
  real x2_pred[N2]; // Incomplete observation covariate retrodiction
  real y2_pred[N2]; // Incomplete observation conditional variate prediction
  
  real f_grid[N_grid]; // Inferred functional behavior along covariate grid
  real y2_pred_grid[N_grid]; // Predictive variate behavior along covariate grid
  
  for (n in 1:N1) {
    x1_pred[n] = gamma_ld_rng(mu1, phi1);
    y1_pred[n] = gamma_ld_rng(lambda * x1[n], psi + 0.001 * gamma * mu1);
  }
  
  for (n in 1:N2) {
    x2_pred[n] = gamma_ld_rng(mu2, phi2);
    y2_pred[n] = gamma_ld_rng(lambda * x2[n], psi + 0.001 * gamma * mu2);
  }
  
  for (n in 1:N_grid) {
    f_grid[n] = lambda * x_grid[n];
    y2_pred_grid[n] = gamma_ld_rng(f_grid[n], psi + 0.001 * gamma * mu2);
  }
}
