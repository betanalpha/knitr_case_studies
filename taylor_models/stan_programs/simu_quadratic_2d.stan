data {
  vector[2] x0;
  real<lower=0> tau_x;
}

transformed data {
  int<lower=0> M = 2;   // Number of covariates
  int<lower=0> N = 250; // Number of observations
  
  vector[2] z0 = [-5, 8]';
  
  real gamma0 = -7;                          // Intercept
  vector[M] gamma1 = [4, -6]';               // Linear slopes
  matrix[M, M] gamma2 = [ [3, -1], [-1, 1]]; // Quadratic slopes
  real sigma = 0.5;                          // Measurement variability
}

generated quantities {  
  matrix[N, M] X; // Covariate design matrix
  real y[N];      // Variates

  for (n in 1:N) {
    X[n, 1] = normal_rng(x0[1], tau_x);
    X[n, 2] = normal_rng(x0[2], tau_x);
    y[n] = normal_rng(  gamma0 
                      + (X[n] - z0') * gamma1
                      + (X[n] - z0') * gamma2 * (X[n] - z0')', sigma);
  }
}
