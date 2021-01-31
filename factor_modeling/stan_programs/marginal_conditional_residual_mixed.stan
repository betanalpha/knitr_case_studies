data {
  int<lower=0> N;      // Number of observations
  vector[N] y;         // Observations
  real<lower=0> sigma; // Measurement variability
  
  // Context information
  int<lower=0> K;
  int<lower=1, upper=K> obs_to_context[N];
  
  int<lower=0, upper=K> K_ncp;                  // Number of noncentered contexts
  int<lower=1, upper=K> context_ncp_idx[K_ncp]; // Index of noncentered contexts
  
  int<lower=0, upper=K> K_cp;                   // Number of centered contexts
  int<lower=1, upper=K> context_cp_idx[K_cp];   // Index of noncentered contexts
  
  // Factor information
  int<lower=0> N_factor_levels;
  int<lower=1, upper=N_factor_levels> context_to_level[K];
  
  int<lower=0, upper=N_factor_levels> L_ncp;                 // Number of noncentered levels
  int<lower=1, upper=N_factor_levels> factor_ncp_idx[L_ncp]; // Index of noncentered levels

  int<lower=0, upper=N_factor_levels> L_cp;                  // Number of centered levels
  int<lower=1, upper=N_factor_levels> factor_cp_idx[L_cp];   // Index of noncentered levels
}

parameters {
  // Hierarchical models for the within-level residuals
  vector[K_ncp] delta_context_ncp;                    // Non-centered residuals
  vector[K_cp] delta_context_cp;                      // Centered residuals
  vector<lower=0>[N_factor_levels] delta_context_tau; // Residual scales
  
  // Hierarchical model for between-level residuals
  vector[L_ncp] delta_factor_ncp; // Non-centered residuals
  vector[L_cp]  delta_factor_cp;  // Centered residuals
  real<lower=0> delta_factor_tau; // Residual scale
  
  // Baseline
  real baseline;
}

transformed parameters {
  vector[K] delta_context;
  vector[N_factor_levels] delta_factor;
  vector[K] theta_context;
  
  delta_context[context_ncp_idx] 
    = delta_context_tau[context_to_level[context_ncp_idx]] .* delta_context_ncp;
  delta_context[context_cp_idx] = delta_context_cp;
  
  delta_factor[factor_ncp_idx] = delta_factor_tau * delta_factor_ncp;
  delta_factor[factor_cp_idx] = delta_factor_cp;
  
  theta_context =  baseline                        // Baseline
                 + delta_factor[context_to_level]  // Between-level variation
                 + delta_context;                  // Within-level variation
}

model {
  delta_context_cp ~ normal(0, delta_context_tau[context_to_level[context_cp_idx]]);
  delta_context_ncp ~ normal(0, 1);
  delta_context_tau ~ normal(0, 3);
    
  delta_factor_cp ~ normal(0, delta_factor_tau);
  delta_factor_ncp ~ normal(0, 1);
  delta_factor_tau ~ normal(0, 3);
  
  baseline ~ normal(0, 5);
  
  y ~ normal(theta_context[obs_to_context], sigma);
}
