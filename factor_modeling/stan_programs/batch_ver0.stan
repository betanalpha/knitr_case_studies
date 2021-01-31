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
  vector[N] log_rho_I; // Observed input signal
  int y[N];            // Observered output counts
  
  // First factor information
  int<lower=0> N_factor1_levels;
  
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
}

transformed parameters {
  vector[N_factor3_levels] log_psi = rep_vector(baseline, N_factor3_levels);
}

model {  
  baseline ~ normal(0, 1);

  y ~ poisson_log(log_psi[obs_to_factor3] + log_rho_I);
}

generated quantities {
  int y_post_pred[N] = poisson_log_rng(log_psi[obs_to_factor3] + log_rho_I);
  real F_post_pred = F_stat(to_vector(y_post_pred) ./ exp(log_rho_I), 
                            N, obs_to_factor3, N_factor3_levels);
}
