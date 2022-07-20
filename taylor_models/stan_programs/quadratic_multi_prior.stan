data {
  int<lower=0> M; // Number of covariates
  int<lower=0> N; // Number of observations
  
  vector[M] x0;   // Covariate baselines
  matrix[N, M] X; // Covariate design matrix
}

transformed data {
  matrix[N, M] deltaX;
  for (n in 1:N) {
    deltaX[n,] = X[n] - x0';
  }
}

parameters {
  real alpha;                      // Intercept
  vector[M] beta1;                 // Linear slopes
  vector[M] beta2_d;               // On-diagonal quadratic slopes
  vector[M * (M - 1) / 2] beta2_o; // Off-diagonal quadratic slopes
  real<lower=0> sigma;             // Measurement Variability
}

model {
  // Prior model
  alpha ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta2_d ~ normal(0, 2);
  beta2_o ~ normal(0, 1);
  sigma ~ normal(0, 5);
}

generated quantities {
  vector[N] f = alpha + deltaX * beta1;
  {
    matrix[M, M] beta2;
    for (m1 in 1:M) {
      beta2[m1, m1] = beta2_d[m1];
      for (m2 in (m1 + 1):M) {
        int m3 = (2 * M - m1) * (m1 - 1) / 2 + m2 - m1;
        beta2[m1, m2] = beta2_o[m3];
        beta2[m2, m1] = beta2_o[m3];
      }
    }
    
    for (n in 1:N) {
      f[n] = f[n] + deltaX[n] * beta2 * deltaX[n]';
    }
  }
}
