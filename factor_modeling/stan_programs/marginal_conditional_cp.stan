data {
  int<lower=0> N;      // Number of observations
  vector[N] y;         // Observations
  real<lower=0> sigma; // Measurement variability
  
  // Number of individual contexts
  int<lower=0> K;
  int<lower=1, upper=K> obs_to_context[N];
  
  // Factor information
  int<lower=0> N_factor_levels;
  int<lower=1, upper=N_factor_levels> context_to_level[K];
}

parameters {
  // Centered hierarchical models for the conditional heterogeneity within each level
  vector[K] context_theta;                       // Centered individual parameters
  vector<lower=0>[N_factor_levels] residual_tau; // Population scales within each level
  
  // Centered hierarchical model for the marginal level heterogeneity
  vector[N_factor_levels] factor_theta; // Centered individual parameters
  real factor_mu;                       // Population location
  real<lower=0> factor_tau;             // Population scale
}

model {
  context_theta ~ normal(factor_theta[context_to_level],
                         residual_tau[context_to_level]);
  residual_tau ~ normal(0, 3);
    
  factor_theta ~ normal(factor_mu, factor_tau);
  factor_mu ~ normal(0, 5);
  factor_tau ~ normal(0, 3);
  
  y ~ normal(context_theta[obs_to_context], sigma);
}
