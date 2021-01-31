functions {
  real F_stat(vector y, int N, int[] indiv_idxs, int N_indiv) {
    real global_mean = 0;
    vector[N_indiv] total_counts = rep_vector(0, N_indiv);
    vector[N_indiv] running_counts = rep_vector(0, N_indiv);
    vector[N_indiv] means = rep_vector(0, N_indiv);
    vector[N_indiv] sum_squares = rep_vector(0, N_indiv);
    
    int N_min_occupunacy = 2; // Minimum occupuancy for included contexts
    real N_occupied = 0;      // Number of contexts with at least two counts
    real between = 0; // Empirical variance of empirical means in occupied contexts
    real within = 0;  // Empirical mean of emperical variances in occupied contexts

    for (n in 1:N) 
      total_counts[indiv_idxs[n]] += 1;
  
    for (n in 1:N) {
      int idx = indiv_idxs[n];
      real delta;
      
      // Skip (nearly) unoccupied contexts
      if (total_counts[idx] < N_min_occupunacy) continue;
      
      delta = y[n] - global_mean;
      global_mean += delta / n;

      running_counts[idx] += 1;
      delta = y[n] - means[idx];
      means[idx] += delta / running_counts[idx];
      sum_squares[idx] += (y[n] - means[idx]) * delta;
    }

    for (l in 1:N_indiv) {
      // Skip (nearly) unoccupied contexts
      if (total_counts[l] < N_min_occupunacy) continue;
      
      N_occupied += 1;
      between += total_counts[l] * square(means[l] - global_mean);
      within += (sum_squares[l] - within) / N_occupied;
    }
    
    // Return zero if not enough minimally occupied contexts
    if (N_occupied < 2) return 0; 
    
    between = between / (N_occupied - 1);
    within /= (N - N_occupied);
    return between / within;
  }
}

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
    vector[N_main_factor_levels[f]] delta_main
      = packed_delta_main[main_factor_start_idx[f]:main_factor_end_idx[f]];
    theta += delta_main[main_factor_level_idx[f]];
  }
  
  packed_delta_main_ncp ~ normal(0, 1);
  packed_delta_main_cp ~ normal(0, delta_main_tau[packed_main_factor_idx[main_cp_idx]]);
  delta_main_tau ~ normal(0, 4);
  
  // First-order factor interaction contributions
  for (f in 1:N_inter1_factors) {
    vector[N_inter1_factor_levels[f]] delta_inter1
      = packed_delta_inter1[inter1_factor_start_idx[f]:inter1_factor_end_idx[f]];
    theta += delta_inter1[inter1_factor_level_idx[f]];
  }
  
  packed_delta_inter1_ncp ~ normal(0, 1);
  packed_delta_inter1_cp ~ normal(0, delta_inter1_tau[packed_inter1_factor_idx[inter1_cp_idx]]);
  delta_inter1_tau ~ normal(0, 2);
  
  // Second-order factor interaction contributions
  for (f in 1:N_inter2_factors) {
    vector[N_inter2_factor_levels[f]] delta_inter2
      = packed_delta_inter2[inter2_factor_start_idx[f]:inter2_factor_end_idx[f]];
    theta += delta_inter2[inter2_factor_level_idx[f]];
  }
  
  packed_delta_inter2_ncp ~ normal(0, 1);
  packed_delta_inter2_cp ~ normal(0, delta_inter2_tau[packed_inter2_factor_idx[inter2_cp_idx]]);
  delta_inter2_tau ~ normal(0, 1);

  // Observational model
  y ~ normal(theta, sigma);
}

generated quantities {
  real y_post_pred[N];
  
  real F_main_post_pred[N_main_factors];
  real F_inter1_post_pred[N_inter1_factors];
  real F_inter2_post_pred[N_inter2_factors];
  
  {
    vector[N] theta =  rep_vector(baseline, N);

    for (f in 1:N_main_factors) {
      vector[N_main_factor_levels[f]] delta_main
        = packed_delta_main[main_factor_start_idx[f]:main_factor_end_idx[f]];
      theta += delta_main[main_factor_level_idx[f]];
    }
    
    for (f in 1:N_inter1_factors) {
      vector[N_inter1_factor_levels[f]] delta_inter1
        = packed_delta_inter1[inter1_factor_start_idx[f]:inter1_factor_end_idx[f]];
      theta += delta_inter1[inter1_factor_level_idx[f]];
    }
    
    for (f in 1:N_inter2_factors) {
      vector[N_inter2_factor_levels[f]] delta_inter2
        = packed_delta_inter2[inter2_factor_start_idx[f]:inter2_factor_end_idx[f]];
      theta += delta_inter2[inter2_factor_level_idx[f]];
    }
    
    y_post_pred = normal_rng(theta, sigma);
  
    for (f in 1:N_main_factors) {
      F_main_post_pred[f] = F_stat(to_vector(y_post_pred), N, 
                                   main_factor_level_idx[f], 
                                   N_main_factor_levels[f]);
    }
    
    for (f in 1:N_inter1_factors) {
      F_inter1_post_pred[f] = F_stat(to_vector(y_post_pred), N, 
                                     inter1_factor_level_idx[f], 
                                     N_inter1_factor_levels[f]);
    }
    
    for (f in 1:N_inter2_factors) {
      F_inter2_post_pred[f] = F_stat(to_vector(y_post_pred), N, 
                                     inter2_factor_level_idx[f], 
                                     N_inter2_factor_levels[f]);
    }
  }
}
