transformed data {
  int N_factor_configs_latent = 250;
  int N_time_obs = 5;
  int N_latent = N_factor_configs_latent * N_time_obs;
  
  int N1 = 9;
  int N2 = 5;
  int N3 = 7;
  
  real alpha = 2;
  real baseline = -4.0 / 365;
  
  // Location
  real<lower=0> delta1_tau = 1.0 / 365;

  // Maintanence Method
  real<lower=0> delta2_tau = 0.5 / 365;

  // Source Quality
  real<lower=0> delta3_tau = 0.005 / 365;

  // First-order interactions
  int<lower=0> N_factor12_levels = N1 * N2;
  real<lower=0> delta12_tau = 0.001 / 365;

  int<lower=0> N_factor13_levels = N1 * N3;
  real<lower=0> delta13_tau = 0.001 / 365;

  int<lower=0> N_factor23_levels = N2 * N3;
  real<lower=0> delta23_tau = 0.5 / 365;

  vector[N1] factor1_level_probs = rep_vector(1.0 / N1, N1);
  vector[N2] factor2_level_probs = rep_vector(1.0 / N2, N2);
  vector[N3] factor3_level_probs = rep_vector(1.0 / N3, N3);
}

generated quantities {
  int N_factor_configs = N_factor_configs_latent;
  int N = N_latent;
  int N_samples = 25;
  
  vector[N1] delta1 = to_vector(normal_rng(rep_vector(0, N1), delta1_tau));
  vector[N2] delta2 = to_vector(normal_rng(rep_vector(0, N2), delta2_tau));
  vector[N3] delta3 = to_vector(normal_rng(rep_vector(0, N3), delta3_tau));
  
  vector[N1 * N2] delta12 = to_vector(normal_rng(rep_vector(0, N1 * N2), delta12_tau));
  vector[N1 * N3] delta13 = to_vector(normal_rng(rep_vector(0, N1 * N3), delta13_tau));
  vector[N2 * N3] delta23 = to_vector(normal_rng(rep_vector(0, N2 * N3), delta23_tau));
  
  int<lower=0> N_main_factors = 3;
  int<lower=0> N_main_factor_levels[3] = {N1, N2, N3};
  int<lower=1> main_factor_level_idx[3, N_factor_configs_latent];
  
  int<lower=1> factor_config_idx[N_latent];
    
  int y[N_latent];
  real t[N_latent];

  { 
    vector[N_factor_configs] kappa;
    
    int factor1_idx[N_factor_configs];
    int factor2_idx[N_factor_configs];
    int factor3_idx[N_factor_configs];
    
    int factor12_idx[N_factor_configs];
    int factor13_idx[N_factor_configs];
    int factor23_idx[N_factor_configs];

    for (n in 1:N_factor_configs) {
      factor1_idx[n] = categorical_rng(factor1_level_probs);
      factor2_idx[n] = categorical_rng(factor2_level_probs);
      factor3_idx[n] = categorical_rng(factor3_level_probs);
      
      factor12_idx[n] = (factor1_idx[n] - 1) * N2 + factor2_idx[n];
      factor23_idx[n] = (factor2_idx[n] - 1) * N3 + factor3_idx[n];
      factor13_idx[n] = (factor1_idx[n] - 1) * N3 + factor3_idx[n];
    }
      
    main_factor_level_idx[1] = factor1_idx;
    main_factor_level_idx[2] = factor2_idx;
    main_factor_level_idx[3] = factor3_idx;

    kappa =  baseline
           + delta1[factor1_idx]
           + delta2[factor2_idx]
           + delta3[factor3_idx]
           + delta12[factor12_idx]
           + delta13[factor13_idx]
           + delta23[factor23_idx];
            
    for (c in 1:N_factor_configs) {
      for (o in 1:N_time_obs) {
        int n = (c - 1) * N_time_obs + o;
        int quarter = o - 1;
        int delta_day = binomial_rng(6, 0.5) - 3;
        if (quarter == 0)
          t[n] = fabs(delta_day);
        else
          t[n] = 90 * quarter + delta_day;
        
        factor_config_idx[n] = c;
        y[n] = binomial_rng(N_samples, inv_logit(alpha + kappa[c] * t[n]));
      }
    }
  }
}
