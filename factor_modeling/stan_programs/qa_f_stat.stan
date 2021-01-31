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
  int<lower=0> N;                     // Number of observations
  int<lower=0> N_samples;             // Number of items sampled at each observation
  int<lower=0, upper=N_samples> y[N]; // Number of passing items at each observation
  vector[N] t;                        // Observation times (days)
  
  // Days for conditional retrodictive checks
  int<lower=0> N_checks;
  int<lower=0> N_check_obs;
  int<lower=1, upper=N> check_idx[N_checks, N_check_obs]; 
  
  // Total factor configuration information
  int<lower=0> N_factor_configs;
  int<lower=0, upper=N_factor_configs> factor_config_idx[N];
  
  // Main factor contributions logistics
  int<lower=0> N_main_factors;
  int<lower=0> N_main_factor_levels[N_main_factors];
  int<lower=1> main_factor_level_idx[N_main_factors, N_factor_configs];
  
  // First-order factor interaction contributions logistics
  int<lower=0> N_inter1_factors;
  int<lower=0> N_inter1_factor_levels[N_inter1_factors];
  int<lower=1> inter1_factor_level_idx[N_inter1_factors, N_factor_configs];

  // Second-order factor interaction contributions logistics
  int<lower=0> N_inter2_factors;
  int<lower=0> N_inter2_factor_levels[N_inter2_factors];
  int<lower=1> inter2_factor_level_idx[N_inter2_factors, N_factor_configs];
}

generated quantities {
  real F_main_obs[N_main_factors, N_checks];
  real F_inter1_obs[N_inter1_factors, N_checks];
  real F_inter2_obs[N_inter2_factors, N_checks];
  
  for (c in 1:N_checks) {
    int check_to_factor[N_check_obs] = factor_config_idx[check_idx[c]];
    
    for (f in 1:N_main_factors)
      F_main_obs[f, c] = F_stat(to_vector(y[check_idx[c]]), N_check_obs, 
                                main_factor_level_idx[f][check_to_factor], 
                                N_main_factor_levels[f]);
      
    for (f in 1:N_inter1_factors)
      F_inter1_obs[f, c] = F_stat(to_vector(y[check_idx[c]]), N_check_obs,
                                  inter1_factor_level_idx[f][check_to_factor], 
                                  N_inter1_factor_levels[f]);
    for (f in 1:N_inter2_factors)
      F_inter2_obs[f, c] = F_stat(to_vector(y[check_idx[c]]), N_check_obs, 
                                  inter2_factor_level_idx[f][check_to_factor], 
                                  N_inter2_factor_levels[f]);  
  }
}
