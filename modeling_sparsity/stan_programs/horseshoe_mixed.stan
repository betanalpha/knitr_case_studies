data {
  int<lower=0> K;                       // Number of groups
  int<lower=0> N;                       // Number of observations
  vector[N] y;                          // Observations
  int<lower=1, upper=K> context_idx[N]; // Context assignments
  
  int<lower=0, upper=K> K_cp;
  int<lower=0, upper=K> cp_idx[K_cp];
  
  int<lower=0, upper=K> K_ncp;
  int<lower=0, upper=K> ncp_idx[K_ncp];
  
  real<lower=0> sigma;                  // Measurement variability
}

parameters {
  vector[K_cp] theta_cp;
  vector[K_ncp] theta_ncp;
  
  vector<lower=0>[K] lambda;
  real<lower=0> tau;
}

transformed parameters {
  vector[K] theta;
  theta[cp_idx] = theta_cp;
  theta[ncp_idx] = theta_ncp .* lambda[ncp_idx] * tau;
}

model {
  // Horseshoe prior model
  theta_cp ~ normal(0, lambda[cp_idx] * tau);
  theta_ncp ~ normal(0, 1);
  
  lambda ~ cauchy(0, 1);
  tau ~ normal(0, 0.1);

  // Observational model
  y ~ normal(theta[context_idx], sigma);
}
