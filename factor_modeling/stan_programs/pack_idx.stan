data {
  // Main factor indices
  int<lower=0> N_main_factors;
  int<lower=0> N_main_factor_levels[N_main_factors];
  
  // First-order interaction indices
  int<lower=0> N_inter1_factors;
  int<lower=0> N_inter1_factor_levels[N_inter1_factors];
  
  // Second-order factor indices
  int<lower=0> N_inter2_factors;
  int<lower=0> N_inter2_factor_levels[N_inter2_factors];
}

transformed data {
  int N_all_main = sum(N_main_factor_levels);
  int N_all_inter1 = sum(N_inter1_factor_levels);
  int N_all_inter2 = sum(N_inter2_factor_levels);
}

generated quantities {
  // Flattened indices
  int<lower=0> N_all_main_factor_levels = N_all_main;
  int<lower=1, upper=N_main_factors> packed_main_factor_idx[N_all_main];
  int<lower=1, upper=N_all_main> main_factor_start_idx[N_main_factors];
  int<lower=1, upper=N_all_main> main_factor_end_idx[N_main_factors];
  
  int<lower=0> N_all_inter1_factor_levels = N_all_inter1;
  int<lower=1, upper=N_inter1_factors> packed_inter1_factor_idx[N_all_inter1];
  int<lower=1, upper=N_all_inter1> inter1_factor_start_idx[N_inter1_factors];
  int<lower=1, upper=N_all_inter1> inter1_factor_end_idx[N_inter1_factors];
  
  int<lower=0> N_all_inter2_factor_levels = N_all_inter2;
  int<lower=1, upper=N_inter2_factors> packed_inter2_factor_idx[N_all_inter2];
  int<lower=1, upper=N_all_inter2> inter2_factor_start_idx[N_inter2_factors];
  int<lower=1, upper=N_all_inter2> inter2_factor_end_idx[N_inter2_factors];

  {
    int cumsum = 0;
    for (f in 1:N_main_factors) {
      main_factor_start_idx[f] = cumsum + 1;
      main_factor_end_idx[f] = cumsum + N_main_factor_levels[f];
      packed_main_factor_idx[main_factor_start_idx[f]:main_factor_end_idx[f]] 
        = rep_array(f, N_main_factor_levels[f]);
      cumsum += N_main_factor_levels[f];
    }
    
    cumsum = 0;
    for (f in 1:N_inter1_factors) {
      inter1_factor_start_idx[f] = cumsum + 1;
      inter1_factor_end_idx[f] = cumsum + N_inter1_factor_levels[f];
      packed_inter1_factor_idx[inter1_factor_start_idx[f]:inter1_factor_end_idx[f]] 
        = rep_array(f, N_inter1_factor_levels[f]);
      cumsum += N_inter1_factor_levels[f];
    }
    
    cumsum = 0;
    for (f in 1:N_inter2_factors) {
      inter2_factor_start_idx[f] = cumsum + 1;
      inter2_factor_end_idx[f] = cumsum + N_inter2_factor_levels[f];
      packed_inter2_factor_idx[inter2_factor_start_idx[f]:inter2_factor_end_idx[f]] 
        = rep_array(f, N_inter2_factor_levels[f]);
      cumsum += N_inter2_factor_levels[f];
    }
  }
}
