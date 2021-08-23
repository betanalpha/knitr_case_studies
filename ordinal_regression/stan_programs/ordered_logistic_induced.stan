functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}

data {
  int<lower=1> N;             // Number of observations
  int<lower=1> K;             // Number of ordinal categories
  int<lower=1, upper=K> y[N]; // Observed ordinals
}

parameters {
  real gamma;       // Latent effect
  ordered[K - 1] c; // (Internal) cut points
}

model {
  // Prior model
  gamma ~ normal(0, 1);
  c ~ induced_dirichlet(rep_vector(1, K), 0);

  // Observational model
  y ~ ordered_logistic(rep_vector(gamma, N), c);
}

generated quantities {
  int<lower=1, upper=K> y_ppc[N];
  for (n in 1:N)
    y_ppc[n] = ordered_logistic_rng(gamma, c);
}
