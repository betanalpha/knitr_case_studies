transformed data {
  int<lower=0> M = 3;         // Number of covariates
  int<lower=0> N = 1000;      // Number of observations
  
  vector[M] x0 = [0, 2, -1]'; // Covariate baseline
  vector[M] z0 = [-1, 3, 0]'; // Latent functional behavior baseline
  real gamma0 = 4.75;                     // True intercept
  vector[M] gamma1 = [-0.4, 0.375, 0.5]'; // True slopes
  real phi = 5;                           // True negative binomial concentration
}

generated quantities {
  matrix[N, M] X; // Covariate design matrix  
  real y[N];      // Variates
  real log_mu[N]; 

  for (n in 1:N) {
    X[n, 1] = normal_rng(x0[1], 0.25);
    X[n, 2] = normal_rng(  x0[2] 
                         - 3 * (X[n, 1] - x0[1])
                         + 1 * pow(X[n, 1] - x0[1], 3), 0.3);
    
    X[n, 3] =     x0[3] 
              + ( (X[n, 1] - x0[1]) - 10 * fabs(normal_rng(0, 2.25)) ) / sqrt(1 + 100);
    
    /*
    if ((X[n, 1] - x0[1]) * (X[n, 2] - x0[2]) > 0)
      X[n, 3] 
        = normal_rng(x0[3] + 0.41 * sqrt((X[n, 1] - x0[1]) * (X[n, 2] - x0[2])), 
                     0.5);
    else
      X[n, 3] 
        = normal_rng(x0[3] - 0.41 * sqrt(-(X[n, 1] - x0[1]) * (X[n, 2] - x0[2])), 
                     0.4);
    */
      
    log_mu[n] = gamma0 + (X[n] - z0') * gamma1;
    y[n] = neg_binomial_2_log_rng(gamma0 + (X[n] - z0') * gamma1, phi);
  }
}
