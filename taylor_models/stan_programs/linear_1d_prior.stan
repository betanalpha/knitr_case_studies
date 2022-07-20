data {
  real x0; // Baseline input
  
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
}

generated quantities {
  real f_grid[N_grid];
  for (n in 1:N_grid) {
    f_grid[n] = alpha + beta * (x_grid[n] - x0);
  }
}
