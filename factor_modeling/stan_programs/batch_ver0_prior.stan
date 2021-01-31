data {
  int<lower=0> N;      // Number of observations
  
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
}
