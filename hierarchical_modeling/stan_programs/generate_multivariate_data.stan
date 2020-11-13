data {
  // Number of observations
  int<lower=1> N;      
  
  // Number of individuals in hierarchy
  int<lower=0> K;
  
  // Individual from which each observation is generated
  int<lower=1, upper=K> indiv_idx[N];
  
   // Number of components in each observation
  int<lower=1> I;     
  
  // Measurement variability
  real<lower=0> sigma; 
}

transformed data {
  // Population locations
  vector[I] mu = [-2.191979, 1.782117, 2.349158]';
  
  // Population scales
  vector<lower=0>[I] tau = [1.0625921, 0.7013441, 0.2720867]';
  
  // Population correlations
  matrix[I, I] L = [ [1.0000000, 0.0000000, 0.0000000],
                     [0.1885095, 0.9820714, 0.0000000],
                     [0.4708125, 0.2414451, 0.8485516] ];
                     
  vector[I] theta[K];
  for (k in 1:K)
    theta[k] = multi_normal_cholesky_rng(mu, diag_pre_multiply(tau, L));
}

generated quantities {
  real y[N, I]; // Simulated observations
  
  for (n in 1:N) {
    // Simulate individual observations
    y[n] = normal_rng(theta[indiv_idx[n]], sigma);
  } 
}
