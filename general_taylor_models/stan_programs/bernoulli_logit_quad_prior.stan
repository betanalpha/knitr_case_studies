data {
  int<lower=0> M; // Number of covariates
  int<lower=0> N; // Number of observations
  
  vector[M] x0;   // Covariate baselines
  matrix[N, M] X; // Covariate design matrix
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
  real alpha;                      // Intercept
  vector[M] beta1;                 // Linear slopes
  vector[M] beta2_d;               // On-diagonal quadratic slopes
  vector[M * (M - 1) / 2] beta2_o; // Off-diagonal quadratic slopes
}

model {
  // Prior model
  alpha ~ normal(0, 1);
  beta1 ~ normal(0, 1);
  beta2_d ~ normal(0, 0.5);
  beta2_o ~ normal(0, 0.5);
}

// Simulate a full observation from the current value of the parameters
generated quantities {
  vector[N] logit_p;
  int y_pred[N];
  {
    vector[M * (M + 3) / 2 + 1] beta
      = append_row(
          append_row(
            append_row(alpha, beta1), 
          beta2_d), 
        beta2_o);
    logit_p = deltaX * beta;
    y_pred = bernoulli_logit_rng(logit_p);
  }
}
