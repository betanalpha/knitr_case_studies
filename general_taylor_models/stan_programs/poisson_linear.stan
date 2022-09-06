data {
  int<lower=0> M; // Number of covariates
  int<lower=0> N; // Number of observations
  
  vector[M] x0;   // Covariate baselines
  matrix[N, M] X; // Covariate design matrix
  
  int y[N];       // Variates
}

transformed data {
  matrix[N, M] deltaX;
  for (n in 1:N) {
    deltaX[n,] = X[n] - x0';
  }
}

parameters {
  real alpha;      // Intercept
  vector[M] beta;  // Slopes
}

model {
  // Prior model
  alpha ~ normal(25, 5);
  beta ~ normal(0, 3);

  // Vectorized observation model
  y ~ poisson(alpha + deltaX * beta);
}

// Simulate a full observation from the current value of the parameters
generated quantities {
  vector[N] mu = alpha + deltaX * beta;
  int y_pred[N] = poisson_rng(mu);
}
