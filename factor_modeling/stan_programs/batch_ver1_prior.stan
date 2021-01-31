data {
  int<lower=0> N;      // Number of observations
  
  // First factor information
  int<lower=0> N_factor1_levels;
  
  int<lower=0, upper=N_factor1_levels> L_factor1_ncp;
  int<lower=1, upper=N_factor1_levels> factor1_ncp_idx[L_factor1_ncp];

  int<lower=0, upper=N_factor1_levels> L_factor1_cp;
  int<lower=1, upper=N_factor1_levels> factor1_cp_idx[L_factor1_cp];
  
  // Second factor information
  int<lower=0> N_factor2_levels;
  
  // Third factor information
  int<lower=0> N_factor3_levels;
  
  // Nesting configuration
  int<lower=1, upper=N_factor3_levels> obs_to_factor3[N];
  int<lower=1, upper=N_factor2_levels> factor3_to_factor2[N_factor3_levels];
  int<lower=1, upper=N_factor1_levels> factor2_to_factor1[N_factor2_levels];
}

parameters {  
  // Baseline
  real baseline;
  
  // Hierarchical model for the first factor residuals
  vector[L_factor1_ncp] delta_factor1_ncp;
  vector[L_factor1_cp] delta_factor1_cp;
  real<lower=0> delta_factor1_tau;
}

transformed parameters {
  vector[N_factor1_levels] delta_factor1;

  vector[N_factor3_levels] log_psi;

  // Reconstruct first factor residuals
  delta_factor1[factor1_ncp_idx] = delta_factor1_tau * delta_factor1_ncp;
  delta_factor1[factor1_cp_idx] = delta_factor1_cp;
  
  log_psi =  baseline
           + delta_factor1[factor2_to_factor1[factor3_to_factor2]];
}

model {  
  baseline ~ normal(0, 1);

  delta_factor1_cp ~ normal(0, delta_factor1_tau);
  delta_factor1_ncp ~ normal(0, 1);
  delta_factor1_tau ~ normal(0, 1);
}
