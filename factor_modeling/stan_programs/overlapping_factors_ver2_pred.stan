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

  // Out of sample prediction data
  int<lower=0> pred_N;
  int<lower=1> pred_main_factor_level_idx[N_main_factors, pred_N];
  int<lower=1> pred_inter1_factor_level_idx[N_inter1_factors, pred_N];
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
}

transformed parameters {
  vector[N_all_main_factor_levels] packed_delta_main;
  vector[N_all_inter1_factor_levels] packed_delta_inter1;
    
  packed_delta_main[main_ncp_idx] 
    = delta_main_tau[packed_main_factor_idx[main_ncp_idx]] .* packed_delta_main_ncp;
  packed_delta_main[main_cp_idx] = packed_delta_main_cp;
  
  packed_delta_inter1[inter1_ncp_idx] 
    = delta_inter1_tau[packed_inter1_factor_idx[inter1_ncp_idx]] .* packed_delta_inter1_ncp;
  packed_delta_inter1[inter1_cp_idx] = packed_delta_inter1_cp;
}

model {
  // Baseline contribution
  vector[N] theta =  rep_vector(baseline, N);
  
  baseline ~ normal(0, 5);
  
  // Main factor contributions
  for (f in 1:N_main_factors) {
    vector[N_main_factor_levels[f]] delta_main
      = packed_delta_main[main_factor_start_idx[f]:main_factor_end_idx[f]];
    theta += delta_main[main_factor_level_idx[f]];
  }
  
  packed_delta_main_ncp ~ normal(0, 1);
  packed_delta_main_cp ~ normal(0, delta_main_tau[packed_main_factor_idx[main_cp_idx]]);
  delta_main_tau ~ normal(0, 3);
  
  // First-order factor interaction contributions
  for (f in 1:N_inter1_factors) {
    vector[N_inter1_factor_levels[f]] delta_inter1
      = packed_delta_inter1[inter1_factor_start_idx[f]:inter1_factor_end_idx[f]];
    theta += delta_inter1[inter1_factor_level_idx[f]];
  }
  
  packed_delta_inter1_ncp ~ normal(0, 1);
  packed_delta_inter1_cp ~ normal(0, delta_inter1_tau[packed_inter1_factor_idx[inter1_cp_idx]]);
  delta_inter1_tau ~ normal(0, 2);

  // Observational model
  y ~ normal(theta, sigma);
}

generated quantities {
  real y_post_pred[N];
  real y_post_out_pred[pred_N];
  {
    vector[N] theta =  rep_vector(baseline, N);
    vector[pred_N] pred_theta =  rep_vector(baseline, pred_N);

    for (f in 1:N_main_factors) {
      vector[N_main_factor_levels[f]] delta_main
        = packed_delta_main[main_factor_start_idx[f]:main_factor_end_idx[f]];
      theta += delta_main[main_factor_level_idx[f]];
      pred_theta += delta_main[pred_main_factor_level_idx[f]];
    }
    
    for (f in 1:N_inter1_factors) {
      vector[N_inter1_factor_levels[f]] delta_inter1
        = packed_delta_inter1[inter1_factor_start_idx[f]:inter1_factor_end_idx[f]];
      theta += delta_inter1[inter1_factor_level_idx[f]];
      pred_theta += delta_inter1[pred_inter1_factor_level_idx[f]];
    }
    
    y_post_pred = normal_rng(theta, sigma);
    y_post_out_pred = normal_rng(pred_theta, sigma);
  }
}
