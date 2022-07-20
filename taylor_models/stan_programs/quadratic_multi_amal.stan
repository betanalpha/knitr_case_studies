data {
  int<lower=0> M; // Number of covariates
  int<lower=0> N; // Number of observations
  
  vector[M] x0;   // Covariate baselines
  matrix[N, M] X; // Covariate design matrix
  
  real y[N];      // Variates
}

transformed data {
  matrix[N, M * (M + 3) / 2 + 1] deltaX;
  for (n in 1:N) {
    deltaX[n, 1] = 1;
    
    for (m1 in 1:M) {
      // Linear pertrubations
      deltaX[n, m1 + 1] = X[n, m1] - x0[m1];
    }
    
    for (m1 in 1:M) {
      // On-diagional quadratic pertrubations
      deltaX[n, M + m1 + 1] 
        = deltaX[n, m1 + 1] * deltaX[n, m1 + 1];
  
      for (m2 in (m1 + 1):M) {
        int m3 = (2 * M - m1) * (m1 - 1) / 2 + m2 - m1;
          
        // Off-diagonal quadratic pertrubations
        // Factor of 2 ensures that beta parameters have the
        // same intepretation as the expanded implementation
        deltaX[n, 2 * M + m3 + 1] 
          = 2 * deltaX[n, m1 + 1] * deltaX[n, m2 + 1];
      }
    }
  }
}

parameters {
  real beta0;                      // Intercept
  vector[M] beta1;                 // Linear slopes
  vector[M] beta2_d;               // On-diagonal quadratic slopes
  vector[M * (M - 1) / 2] beta2_o; // Off-diagonal quadratic slopes
  real<lower=0> sigma;             // Measurement Variability
}

model {
  vector[M * (M + 3) / 2 + 1] beta
    = append_row(
        append_row(
          append_row(beta0, beta1), 
        beta2_d), 
      beta2_o);
  
  // Prior model
  beta0 ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta2_d ~ normal(0, 2);
  beta2_o ~ normal(0, 1);
  sigma ~ normal(0, 5);

  // Vectorized observation model
  y ~ normal(deltaX * beta, sigma);
}

// Simulate a full observation from the current value of the parameters
generated quantities {
  real y_pred[N];
  {
    vector[M * (M + 3) / 2 + 1] beta
      = append_row(
          append_row(
            append_row(beta0, beta1), 
          beta2_d),
        beta2_o);
    y_pred = normal_rng(deltaX * beta, sigma);
  }
}
