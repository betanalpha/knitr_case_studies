data {
  int<lower=0> M; // Number of covariates
  int<lower=0> N; // Number of observations
  
  vector[M] x0;   // Covariate baselines
  matrix[N, M] X; // Covariate design matrix
  
  real y[N];      // Variates
}

transformed data {
  matrix[N, M + 1] deltaX;
  for (n in 1:N) {
    deltaX[n, 1] = 1;
    deltaX[n, 2:(M + 1)] = X[n] - x0';
  }
}

parameters {
  real beta0;          // Intercept
  vector[M] beta1;     // Slopes
  real<lower=0> sigma; // Measurement Variability
}

model {
  vector[M + 1] beta = append_row(beta0, beta1);
  
  // Prior model
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  sigma ~ normal(0, 5);

  // Vectorized observation model
  y ~ normal(deltaX * beta, sigma);
}

// Simulate a full observation from the current value of the parameters
generated quantities {
  real y_pred[N] = normal_rng(deltaX * append_row(beta0, beta1), sigma);
}
