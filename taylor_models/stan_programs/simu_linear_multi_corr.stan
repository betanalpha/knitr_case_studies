data {
  vector[3] x0;
}

transformed data {
  int<lower=0> M = 3;           // Number of covariates
  int<lower=0> N = 100;         // Number of observations
  
  vector[M] z0 = [-1, 3, 0]';
  real gamma0 = 10;               // True intercept
  vector[M] gamma1 = [2, -3, 1]'; // True slopes
  real sigma = 0.5;              // True measurement variability
}

generated quantities {
  matrix[N, M] X; // Covariate design matrix
  real y[N];      // Variates

  for (n in 1:N) {
    X[n, 1] = normal_rng(x0[1], 2);
    X[n, 2] = normal_rng(x0[2] + 3 * sin(X[n, 1] - x0[1]), 1.25);
    X[n, 3] = normal_rng(x0[3] + (X[n, 1] - x0[1]), 1);
    y[n] = normal_rng(gamma0 + (X[n] - z0') * gamma1, sigma);
  }
}
