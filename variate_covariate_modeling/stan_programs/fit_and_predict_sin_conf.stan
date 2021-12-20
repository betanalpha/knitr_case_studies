data {
  int<lower=1> N1; // Number of complete observations
  vector[N1] x1;   // Complete observations
  vector[N1] y1;   // Complete observations
  real<lower=0> alpha1; // Complete observation covariate distribution shape
  real<lower=0> beta1;  // Complete observation covariate distribution rate

  int<lower=1> N2;      // Number of incomplete observations
  vector[N2] x2;        // Incomplete observations
  real<lower=0> alpha2; // Incomplete observation covariate distribution shape
  real<lower=0> beta2;  // Incomplete observation covariate distribution rate
  
  int<lower=1> N_grid; // Number of grid points for quantifying functional behavior
  vector[N_grid] x_grid; // Grid points for quantifying functional behavior
}

parameters { 
  real<lower=0> A;                   // Variate units
  real<lower=0, upper=2 * pi()> phi; // Radians
  real<lower=0> omega1;              // Radians per second
  real<lower=0> omega2;              // Radians per second
  real<lower=0> gamma;               // Seconds
  real<lower=0> sigma;               // Variate units
}

model {
  A ~ normal(0, 5);
  // Implicit uniform prior density for phi
  omega1 ~ normal(0, 5);
  omega2 ~ normal(0, 5);
  gamma ~ normal(0, 1);
  sigma ~ normal(0, 0.25);
  
  x1 ~ gamma(alpha1, beta1 + gamma * omega1);
  y1 ~ normal(A * sin(omega1 * x1 + phi), sigma);
  
  x2 ~ gamma(alpha2, beta2 + gamma * omega2);
}

generated quantities {
  real x1_pred[N1]; // Complete observation covariate retrodiction
  real y1_pred[N1]; // Complete observation conditional variate retrodiction
  real x2_pred[N2]; // Incomplete observation covariate retrodiction
  real y2_pred[N2]; // Incomplete observation conditional variate prediction
  
  real f1_grid[N_grid];      // Inferred functional behavior along covariate grid
  real f2_grid[N_grid];      // Inferred functional behavior along covariate grid
  real y2_pred_grid[N_grid]; // Predictive variate behavior along covariate grid
  
  for (n in 1:N1) {
    x1_pred[n] = gamma_rng(alpha1, beta1 + gamma * omega1);
    y1_pred[n] = normal_rng(A * sin(omega1 * x1[n] + phi), sigma);
  }
  
  for (n in 1:N2) {
    x2_pred[n] = gamma_rng(alpha2, beta2 + gamma * omega2);
    y2_pred[n] = normal_rng(A * sin(omega2 * x2[n] + phi), sigma);
  }
  
  for (n in 1:N_grid) {
    f1_grid[n] = A * sin(omega1 * x_grid[n] + phi);
    f2_grid[n] = A * sin(omega2 * x_grid[n] + phi);
    y2_pred_grid[n] = normal_rng(f2_grid[n], sigma);
  }
}
