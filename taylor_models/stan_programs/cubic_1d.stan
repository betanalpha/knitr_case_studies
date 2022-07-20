data {
  real x0; // Baseline input
  
  int<lower=1> N; // Number of observations
  vector[N] x; // Observed inputs
  vector[N] y; // Observed outputs
  
  int<lower=1> N_grid; // Number of grid points
  real x_grid[N_grid];  // Input grid
}

transformed data {
  vector[N] delta_x;
  vector[N] delta_x2;
  vector[N] delta_x3;
  
  for (n in 1:N) {
    delta_x[n] = x[n] - x0;
    delta_x2[n] = pow(delta_x[n], 2.0);
    delta_x3[n] = pow(delta_x[n], 3.0);
  }
}

parameters { 
  real alpha;  // Intercept
  real beta1;  // First-order slope
  real beta2;  // Second-order slope
  real beta3;  // Third-order slope
  real<lower=0> sigma; // Measurement Variability
}

model {
  alpha ~ normal(0, 3);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  y ~ normal(  alpha 
             + beta1 * delta_x 
             + beta2 * delta_x2
             + beta3 * delta_x3, sigma);
}

generated quantities {
  real y_pred[N];
  real f_grid[N_grid];
  real y_pred_grid[N_grid];
  
  for (n in 1:N)
    y_pred[n] = normal_rng(  alpha 
                           + beta1 * delta_x[n] 
                           + beta2 * delta_x2[n]
                           + beta3 * delta_x3[n], sigma);
    
  for (n in 1:N_grid) {
    f_grid[n] =   alpha 
                + beta1 * (x_grid[n] - x0) 
                + beta2 * pow(x_grid[n] - x0, 2.0)
                + beta3 * pow(x_grid[n] - x0, 3.0);
    y_pred_grid[n] = normal_rng(f_grid[n], sigma);
  }
}
