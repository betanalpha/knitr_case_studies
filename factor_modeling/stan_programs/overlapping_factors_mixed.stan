data {
  int<lower=0> N;      // Number of observations
  vector[N] y;         // Observations
  real<lower=0> sigma; // Measurement variability
  
  // Main factor contributions logistics
  int<lower=0> N_main_factors;
  int<lower=0> N_main_factor_levels[N_main_factors];
  int<lower=1> main_factor_level_idx[N_main_factors, N];
  
  int<lower=0> N_all_main_factor_levels;
  int<lower=1, upper=N_main_factors> packed_main_factor_idx[N_all_main_factor_levels];
  int<lower=1, upper=N_all_main_factor_levels> main_factor_start_idx[N_main_factors];
  int<lower=1, upper=N_all_main_factor_levels> main_factor_end_idx[N_main_factors];
  
  int<lower=0, upper=N_all_main_factor_levels> L_main_ncp;
  int<lower=1, upper=N_all_main_factor_levels> main_ncp_idx[L_main_ncp];
  
  int<lower=0, upper=N_all_main_factor_levels> L_main_cp;
  int<lower=1, upper=N_all_main_factor_levels> main_cp_idx[L_main_cp];
  
  // First-order factor interaction contributions logistics
  int<lower=0> N_inter1_factors;
  int<lower=0> N_inter1_factor_levels[N_inter1_factors];
  int<lower=1> inter1_factor_level_idx[N_inter1_factors, N];
  
  int<lower=0> N_all_inter1_factor_levels;
  int<lower=1, upper=N_inter1_factors> packed_inter1_factor_idx[N_all_inter1_factor_levels];
  int<lower=1, upper=N_all_inter1_factor_levels> inter1_factor_start_idx[N_inter1_factors];
  int<lower=1, upper=N_all_inter1_factor_levels> inter1_factor_end_idx[N_inter1_factors];
  
  int<lower=0, upper=N_all_inter1_factor_levels> L_inter1_ncp;
  int<lower=1, upper=N_all_inter1_factor_levels> inter1_ncp_idx[L_inter1_ncp];
  
  int<lower=0, upper=N_all_inter1_factor_levels> L_inter1_cp;
  int<lower=1, upper=N_all_inter1_factor_levels> inter1_cp_idx[L_inter1_cp];

  // Second-order factor interaction contributions logistics
  int<lower=0> N_inter2_factors;
  int<lower=0> N_inter2_factor_levels[N_inter2_factors];
  int<lower=1> inter2_factor_level_idx[N_inter2_factors, N];
  
  int<lower=0> N_all_inter2_factor_levels;
  int<lower=1, upper=N_inter2_factors> packed_inter2_factor_idx[N_all_inter2_factor_levels];
  int<lower=1, upper=N_all_inter2_factor_levels> inter2_factor_start_idx[N_inter2_factors];
  int<lower=1, upper=N_all_inter2_factor_levels> inter2_factor_end_idx[N_inter2_factors];
  
  int<lower=0, upper=N_all_inter2_factor_levels> L_inter2_ncp;
  int<lower=1, upper=N_all_inter2_factor_levels> inter2_ncp_idx[L_inter2_ncp];
  
  int<lower=0, upper=N_all_inter2_factor_levels> L_inter2_cp;
  int<lower=1, upper=N_all_inter2_factor_levels> inter2_cp_idx[L_inter2_cp];
}

parameters {  
  // Zeroth-order contribution (Baseline)
  real baseline;
  
  // Main factor contributions
  vector[L_main_ncp] packed_delta_main_ncp;
  vector[L_main_cp] packed_delta_main_cp; 
  vector<lower=0>[N_main_factors] delta_main_tau;
  
  // First-order factor interaction contributions
  vector[L_inter1_ncp] packed_delta_inter1_ncp;
  vector[L_inter1_cp] packed_delta_inter1_cp; 
  vector<lower=0>[N_inter1_factors] delta_inter1_tau;
  
  // Second-order factor interaction contributions
  vector[L_inter2_ncp] packed_delta_inter2_ncp;
  vector[L_inter2_cp] packed_delta_inter2_cp; 
  vector<lower=0>[N_inter2_factors] delta_inter2_tau;
}

transformed parameters {
  vector[N_all_main_factor_levels] packed_delta_main;
  vector[N_all_inter1_factor_levels] packed_delta_inter1;
  vector[N_all_inter2_factor_levels] packed_delta_inter2;
  
  packed_delta_main[main_ncp_idx] 
    = delta_main_tau[packed_main_factor_idx[main_ncp_idx]] .* packed_delta_main_ncp;
  packed_delta_main[main_cp_idx] = packed_delta_main_cp;
  
  packed_delta_inter1[inter1_ncp_idx] 
    = delta_inter1_tau[packed_inter1_factor_idx[inter1_ncp_idx]] .* packed_delta_inter1_ncp;
  packed_delta_inter1[inter1_cp_idx] = packed_delta_inter1_cp;
  
  packed_delta_inter2[inter2_ncp_idx] 
    = delta_inter2_tau[packed_inter2_factor_idx[inter2_ncp_idx]] .* packed_delta_inter2_ncp;
  packed_delta_inter2[inter2_cp_idx] = packed_delta_inter2_cp;
}

model {
  // Baseline contribution
  vector[N] theta =  rep_vector(baseline, N);
  baseline ~ normal(0, 8);
  
  // Main factor contributions
  for (f in 1:N_main_factors) {
    // Extract level parameters of fth main factor
    vector[N_main_factor_levels[f]] delta_main
      = packed_delta_main[main_factor_start_idx[f]:main_factor_end_idx[f]];
    
    // Add contribution to total heterogeneity
    theta += delta_main[main_factor_level_idx[f]];
  }
  
  // Main factor hierarchial models
  packed_delta_main_ncp ~ normal(0, 1);
  packed_delta_main_cp ~ normal(0, delta_main_tau[packed_main_factor_idx[main_cp_idx]]);
  delta_main_tau ~ normal(0, 4);

  // First-order factor interaction contributions
  for (f in 1:N_inter1_factors) {
    // Extract level parameters of fth first-order interaction factor
    vector[N_inter1_factor_levels[f]] delta_inter1
      = packed_delta_inter1[inter1_factor_start_idx[f]:inter1_factor_end_idx[f]];
      
    // Add contribution to total heterogeneity
    theta += delta_inter1[inter1_factor_level_idx[f]];
  }
  
  // First-order interaction factor hierarchial models
  packed_delta_inter1_ncp ~ normal(0, 1);
  packed_delta_inter1_cp ~ normal(0, delta_inter1_tau[packed_inter1_factor_idx[inter1_cp_idx]]);
  delta_inter1_tau ~ normal(0, 2);
  
  // Second-order factor interaction contributions
  for (f in 1:N_inter2_factors) {
    // Extract level parameters of fth second-order interaction factor
    vector[N_inter2_factor_levels[f]] delta_inter2
      = packed_delta_inter2[inter2_factor_start_idx[f]:inter2_factor_end_idx[f]];
      
    // Add contribution to total heterogeneity
    theta += delta_inter2[inter2_factor_level_idx[f]];
  }
  
  // Second-order interaction factor hierarchial models
  packed_delta_inter2_ncp ~ normal(0, 1);
  packed_delta_inter2_cp ~ normal(0, delta_inter2_tau[packed_inter2_factor_idx[inter2_cp_idx]]);
  delta_inter2_tau ~ normal(0, 1);

  // Observational model
  y ~ normal(theta, sigma);
}
