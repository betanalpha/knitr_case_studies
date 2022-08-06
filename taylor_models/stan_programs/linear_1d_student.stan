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
  real<lower=0> nu_inv; // Inverse degrees of freedom
  real<lower=0> sigma;  // Measurement Variability
}

model {
  alpha ~ normal(0, 3);
  beta ~ normal(0, 1);
  nu_inv ~ normal(0, 0.5);
  sigma ~ normal(0, 1);
  
  y ~ student_t(inv(nu_inv), alpha + beta * (x - x0), sigma);
}

generated quantities {
  real y_pred[N];
  real mu_grid[N_grid];
  real y_pred_grid[N_grid];
  
  for (n in 1:N)
    y_pred[n] = student_t_rng(inv(nu_inv), alpha + beta * (x[n] - x0), sigma);
    
  for (n in 1:N_grid) {
    mu_grid[n] = alpha + beta * (x_grid[n] - x0);
    y_pred_grid[n] = student_t_rng(inv(nu_inv), mu_grid[n], sigma);
  }
}
