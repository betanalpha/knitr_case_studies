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
  
  int<lower=0, upper=N_factor2_levels> L_factor2_ncp;
  int<lower=1, upper=N_factor2_levels> factor2_ncp_idx[L_factor2_ncp];

  int<lower=0, upper=N_factor2_levels> L_factor2_cp;
  int<lower=1, upper=N_factor2_levels> factor2_cp_idx[L_factor2_cp];
  
  // Third factor information
  int<lower=0> N_factor3_levels;
  
  int<lower=0, upper=N_factor3_levels> L_factor3_ncp;
  int<lower=1, upper=N_factor3_levels> factor3_ncp_idx[L_factor3_ncp];

  int<lower=0, upper=N_factor3_levels> L_factor3_cp;
  int<lower=1, upper=N_factor3_levels> factor3_cp_idx[L_factor3_cp];
  
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
  
  // Hierarchical models for the second factor residuals
  vector[L_factor2_ncp] delta_factor2_ncp;
  vector[L_factor2_cp]  delta_factor2_cp;
  vector<lower=0>[N_factor1_levels] delta_factor2_tau;
  
  // Hierarchical models for the third factor residuals
  vector[L_factor3_ncp] delta_factor3_ncp;
  vector[L_factor3_cp]  delta_factor3_cp;
  vector<lower=0>[N_factor2_levels] delta_factor3_tau;
}

transformed parameters {
  vector[N_factor1_levels] delta_factor1;
  vector[N_factor2_levels] delta_factor2;
  vector[N_factor3_levels] delta_factor3;

  vector[N_factor3_levels] log_psi;

  // Reconstruct first factor residuals
  delta_factor1[factor1_ncp_idx] = delta_factor1_tau * delta_factor1_ncp;
  delta_factor1[factor1_cp_idx] = delta_factor1_cp;

  // Reconstruct second factor residuals
  delta_factor2[factor2_ncp_idx] 
    =   delta_factor2_tau[factor2_to_factor1[factor2_ncp_idx]] .* delta_factor2_ncp;
  delta_factor2[factor2_cp_idx] = delta_factor2_cp;

  // Reconstruct third factor residuals
  delta_factor3[factor3_ncp_idx] 
    =   delta_factor3_tau[factor3_to_factor2[factor3_ncp_idx]] .* delta_factor3_ncp;
  delta_factor3[factor3_cp_idx] = delta_factor3_cp;
  
  log_psi =  baseline
           + delta_factor1[factor2_to_factor1[factor3_to_factor2]]
           + delta_factor2[factor3_to_factor2]
           + delta_factor3;
}

model {  
  baseline ~ normal(0, 1);

  delta_factor1_cp ~ normal(0, delta_factor1_tau);
  delta_factor1_ncp ~ normal(0, 1);
  delta_factor1_tau ~ normal(0, 1);
    
  delta_factor2_cp ~ normal(0, delta_factor2_tau[factor2_to_factor1[factor2_cp_idx]]);
  delta_factor2_ncp ~ normal(0, 1);
  delta_factor2_tau ~ normal(0, 0.5);
  
  delta_factor3_cp ~ normal(0, delta_factor3_tau[factor3_to_factor2[factor3_cp_idx]]);
  delta_factor3_ncp ~ normal(0, 1);
  delta_factor3_tau ~ normal(0, 0.25);
}
