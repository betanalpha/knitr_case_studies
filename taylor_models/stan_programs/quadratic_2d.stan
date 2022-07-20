data {
  int<lower=0> N; // Number of observations
  
  vector[2] x0;   // Covariate baselines
  matrix[N, 2] X; // Covariate design matrix
  
  real y[N];      // Variates
}

transformed data {
  matrix[N, 2] deltaX;
  for (n in 1:N) {
    deltaX[n,] = X[n] - x0';
  }
}

parameters {
  real alpha;           // Intercept
  vector[2] beta1;      // Linear slopes
  vector[2] beta2_d;    // On-diagonal quadratic slopes
  real beta2_o;         // Off-diagonal quadratic slopes
  real<lower=0> sigma;  // Measurement Variability
}

transformed parameters {
  vector[N] f = alpha + deltaX * beta1;
  {
    matrix[2, 2] beta2 = [ [beta2_d[1], beta2_o], [beta2_o, beta2_d[2]] ];
    for (n in 1:N) {
      f[n] = f[n] + deltaX[n] * beta2 * deltaX[n]';
    }
  }
}

model {
  // Prior model
  alpha ~ normal(0, 10);
  beta1 ~ normal(0, 10);
  beta2_d ~ normal(0, 2);
  beta2_o ~ normal(0, 1);
  sigma ~ normal(0, 5);

  // Vectorized observation model
  y ~ normal(f, sigma);
}

// Simulate a full observation from the current value of the parameters
generated quantities {
  real y_pred[N] = normal_rng(f, sigma);
  real det_beta2;
  vector[2] x_star;
  {
    matrix[2, 2] beta2 = [ [beta2_d[1], beta2_o], [beta2_o, beta2_d[2]] ];
    det_beta2 = determinant(beta2);
    if (det_beta2 != 0)
      x_star = -0.5 * mdivide_left(beta2, beta1) + x0;
  }
}
