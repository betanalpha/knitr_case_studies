transformed data {
  int<lower=0> M = 3;         // Number of covariates
  int<lower=0> N = 1000;      // Number of observations
  
  vector[M] x0 = [-1, 0, 1]'; // Covariate baseline
  vector[M] z0 = [-3, 1, 2]'; // Latent functional behavior baseline
  real gamma0 = -2.6;                      // True intercept
  vector[M] gamma1 = [0.2, -2.0, 0.33]';   // True slopes
  matrix[M, M] gamma2 = [ [+0.40, -0.05, -0.20],
                          [-0.05, -1.00, -0.05],
                          [-0.20, -0.05, +0.50] ];
}

generated quantities {
  matrix[N, M] X; // Covariate design matrix
  real y[N];      // Variates

  for (n in 1:N) {
    real x2 = -5;
    while (x2 < x0[2] - 4 || x2 > x0[2] + 4)
      x2 = normal_rng(x0[2], 2);
    
    X[n, 2] = x2;
    X[n, 1] = normal_rng(x0[1] + 1.0 * cos(1.5 * (X[n, 2] - x0[2])), 0.3);
    X[n, 3] = normal_rng(x0[3] + 0.76 * (X[n, 1] - x0[1]), 0.5);

    y[n] = bernoulli_logit_rng(  gamma0 
                               + (X[n] - z0') * gamma1
                               + (X[n] - z0') * gamma2 * (X[n] - z0')');
  }
}
