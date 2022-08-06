data {
  real x0; // Baseline input
  
  int<lower=1> N_grid; // Number of grid points
  real x_grid[N_grid];  // Input grid
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
}

generated quantities {
  real mu_grid[N_grid];
  for (n in 1:N_grid) {
    mu_grid[n] =   alpha 
                 + beta1 * (x_grid[n] - x0) 
                 + beta2 * pow(x_grid[n] - x0, 2.0)
                 + beta3 * pow(x_grid[n] - x0, 3.0);
  }
}
