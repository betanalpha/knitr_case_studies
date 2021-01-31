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
  // Centered hierarchical models for the within-level residuals
  vector[K] delta_context;                          // Centered residuals
  real<lower=0> delta_context_tau[N_factor_levels]; // Residual population scales
  
  // Centered hierarchical model for between level residuals
  vector[N_factor_levels] delta_factor; // Centered residuals
  real<lower=0> delta_factor_tau;       // Residual population scales
  
  // Baseline
  real baseline;
}

transformed parameters {
  vector[K] theta_context =  baseline                        // Baseline
                           + delta_factor[context_to_level]  // Between level variation
                           + delta_context;                  // Within level variation
}

model {
  delta_context ~ normal(0, delta_context_tau[context_to_level]);
  delta_context_tau ~ normal(0, 3);
    
  delta_factor ~ normal(0, delta_factor_tau);
  delta_factor_tau ~ normal(0, 3);
  
  baseline ~ normal(0, 5);
  
  y ~ normal(theta_context[obs_to_context], sigma);
}
