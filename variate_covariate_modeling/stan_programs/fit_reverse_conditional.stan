data {
  int<lower=1> N;
  vector[N] x1; // Complete observation
  vector[N] y1; // Complete observation
  vector[N] y2; // Incomplete observation
  
  int<lower=1> N_grid; // Number of grid points for quantifying functional behavior
  vector[N_grid] y_grid; // Grid points for quantifying functional behavior
}

parameters { 
  real mu1;
  real<lower=0> tau1;
  real<lower=0> sigma;
}

model {
  real sigma_sq = square(sigma);
  real tau1_sq = square(tau1);
  real sum_sq_1 = tau1_sq + sigma_sq;
  
  mu1 ~ normal(0, 3);
  tau1 ~ normal(0, 1);
  sigma ~ normal(0, 1);
  
  x1 ~ normal((y1 * tau1_sq + mu1 * sigma_sq) / sum_sq_1, 
              sqrt(tau1_sq * sigma_sq / sum_sq_1));
}

generated quantities {
  // Predict completion of second observation using
  // conditional variate model from first observation
  vector[N_grid] f_grid;  // Inferred functional behavior along covariate grid
  real x2_pred_grid[N_grid]; // Predictive variate behavior along covariate grid
  
  real x2_pred[N];
  {
    real sigma_sq = square(sigma);
    real tau1_sq = square(tau1);
    real sum_sq_1 = tau1_sq + sigma_sq;
    
    f_grid = (y_grid * tau1_sq + mu1 * sigma_sq) / sum_sq_1;
    
    x2_pred_grid = normal_rng((y_grid * tau1_sq + mu1 * sigma_sq) / sum_sq_1, 
                              sqrt(tau1_sq * sigma_sq / sum_sq_1));
    
    x2_pred = normal_rng((y2 * tau1_sq + mu1 * sigma_sq) / sum_sq_1, 
                         sqrt(tau1_sq * sigma_sq / sum_sq_1));
  }
}
