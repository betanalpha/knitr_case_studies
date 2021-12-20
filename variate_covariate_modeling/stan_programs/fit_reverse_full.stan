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
  
  real mu2;
  real<lower=0> tau2;
  
  real<lower=0> sigma;
}

model {
  real sigma_sq = square(sigma);
  
  real tau1_sq = square(tau1);
  real sum_sq_1 = tau1_sq + sigma_sq;
  
  real tau2_sq = square(tau2);
  real sum_sq_2 = tau2_sq + sigma_sq;
  
  mu1 ~ normal(0, 3);
  tau1 ~ normal(0, 1);
  
  mu2 ~ normal(0, 3);
  tau2 ~ normal(0, 1);
  
  sigma ~ normal(0, 1);
  
  y1 ~ normal(mu1, sqrt(sum_sq_1));
  x1 ~ normal((y1 * tau1_sq + mu1 * sigma_sq) / sum_sq_1, 
              sqrt(tau1_sq * sigma_sq / sum_sq_1));
              
  y2 ~ normal(mu2, sqrt(sum_sq_2));
}

generated quantities {
  // Predict completion of second observation
  vector[N_grid] f1_grid;    // Inferred functional behavior along covariate grid
  vector[N_grid] f2_grid;    // Inferred functional behavior along covariate grid
  real x2_pred_grid[N_grid]; // Predictive variate behavior along covariate grid
  
  real x2_pred[N];
  
  {
    real sigma_sq = square(sigma);
    
    real tau1_sq = square(tau1);
    real sum_sq_1 = tau1_sq + sigma_sq;
    
    real tau2_sq = square(tau2);
    real sum_sq_2 = tau2_sq + sigma_sq;
    
    f1_grid = (y_grid * tau1_sq + mu1 * sigma_sq) / sum_sq_1;
    f2_grid = (y_grid * tau2_sq + mu2 * sigma_sq) / sum_sq_2;
    
    x2_pred_grid = normal_rng((y_grid * tau2_sq + mu2 * sigma_sq) / sum_sq_2, 
                              sqrt(tau2_sq * sigma_sq / sum_sq_2));
    
    x2_pred = normal_rng((y2 * tau2_sq + mu2 * sigma_sq) / sum_sq_2, 
                         sqrt(tau2_sq * sigma_sq / sum_sq_2));
  }
}
