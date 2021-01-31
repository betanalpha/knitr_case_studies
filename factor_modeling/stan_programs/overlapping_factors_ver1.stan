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
}

parameters {  
  // Zeroth-order contribution (Baseline)
  real baseline;
  
  // Main factor contributions
  vector[L_main_ncp] packed_delta_main_ncp;
  vector[L_main_cp] packed_delta_main_cp; 
  vector<lower=0>[N_main_factors] delta_main_tau;
}

transformed parameters {
  vector[N_all_main_factor_levels] packed_delta_main;
  
  packed_delta_main[main_ncp_idx] 
    = delta_main_tau[packed_main_factor_idx[main_ncp_idx]] .* packed_delta_main_ncp;
  packed_delta_main[main_cp_idx] = packed_delta_main_cp;
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
  
  // Observational model
  y ~ normal(theta, sigma);
}

generated quantities {
  real y_post_pred[N];
  
  real F_main_post_pred[N_main_factors];
  real F_inter1_post_pred[N_inter1_factors];
  
  {
    vector[N] theta =  rep_vector(baseline, N);

    for (f in 1:N_main_factors) {
      vector[N_main_factor_levels[f]] delta_main
        = packed_delta_main[main_factor_start_idx[f]:main_factor_end_idx[f]];
      theta += delta_main[main_factor_level_idx[f]];
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
  }
}
