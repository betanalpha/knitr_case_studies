data {
  real x0; // Baseline input
  
  int<lower=1> N; // Number of observations
  vector[N] x; // Observed inputs
  vector[N] y; // Observed outputs
  
  int<lower=1> N_grid; // Number of grid points
  real x_grid[N_grid];  // Input grid
}

parameters { 
  real alpha; // Intercept
  real beta;  // Slope
  real<lower=0> sigma; // Measurement Variability
}

model {
  alpha ~ normal(0, 3);
  beta ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  y ~ normal(alpha + beta * (x - x0), sigma);
}

generated quantities {
  real y_pred[N];
  real mu_grid[N_grid];
  real y_pred_grid[N_grid];
  
  for (n in 1:N)
    y_pred[n] = normal_rng(alpha + beta * (x[n] - x0), sigma);
    
  for (n in 1:N_grid) {
    mu_grid[n] = alpha + beta * (x_grid[n] - x0);
    y_pred_grid[n] = normal_rng(mu_grid[n], sigma);
  }
}
