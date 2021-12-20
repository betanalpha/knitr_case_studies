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
  real<lower=0> A;                   // Variate units
  real<lower=0, upper=2 * pi()> phi; // Radians
  real<lower=0> omega;               // Radians per second
  real<lower=0> sigma;               // Variate units
}

model {
  A ~ normal(0, 5);
  // Implicit uniform prior density for phi
  omega ~ normal(0, 5);
  sigma ~ normal(0, 1);

  y1 ~ normal(A * sin(omega * x1 + phi), sigma);
}

generated quantities {
  real y1_pred[N1]; // Complete observation conditional variate retrodiction
  real y2_pred[N2]; // Incomplete observation conditional variate prediction
  
  real f_grid[N_grid];       // Inferred functional behavior along covariate grid
  real y2_pred_grid[N_grid]; // Predictive variate behavior along covariate grid

  for (n in 1:N1) {
    y1_pred[n] = normal_rng(A * sin(omega * x1[n] + phi), sigma);
  }
  
  for (n in 1:N2) {
    y2_pred[n] = normal_rng(A * sin(omega * x2[n] + phi), sigma);
  }
  
  for (n in 1:N_grid) {
    f_grid[n] = A * sin(omega * x_grid[n] + phi);
    y2_pred_grid[n] = normal_rng(f_grid[n], sigma);
  }
}
