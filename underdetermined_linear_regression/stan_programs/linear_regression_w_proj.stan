data {
  int<lower=1> N; // Number of observations
  int<lower=1> M; // Number of covariates
  matrix[N, M] X; // Design matrix
  vector[N] y;    // Observations
}

transformed data {
  // Buffered design matrix to incorporate intercept
  matrix[N, M + 1] X_buff = append_col(X, rep_vector(1.0, N));
  matrix[M + 1, M + 1] Z = X_buff' * X_buff;
  
  // Vector from origin to underdetermined hyperplane
  vector[M + 1] R = Z \ (X_buff' * y);
  
  // Normalization of projection operator onto underdetermined hyperplane
  matrix[N, N] Pnorm = inverse(X_buff * (Z \ X_buff') );
}

parameters {
  vector[M] beta;
  real alpha;
  real<lower=0> sigma;
}

model {
  // Weakly informative containment priors
  beta ~ normal(0, 10);
  alpha ~ normal(0, 2);
  sigma ~ normal(0, 2);

  // Observational model
  // Equivalent to y ~ normal(X * beta + alpha, sigma);
  y ~ normal(X_buff * append_row(beta, alpha), sigma);
}

generated quantities {
  // Projection of parameters orthogonal to underdetermined hyperplane
  vector[M + 1] gamma_perp;

  // Projection of parameters along underdetermined hyperplane
  vector[M + 1] gamma_par;
  
  // Distance from parameters to underdetermined hyperplane
  real delta;
  
  {
    vector[N] z = Pnorm * X_buff * (append_row(beta, alpha) - R);
    gamma_perp = Z \ (X_buff' * z);
    gamma_par = (append_row(beta, alpha) - R) - gamma_perp;
    delta = sqrt(dot_self(gamma_perp));
  }
}
