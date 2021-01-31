transformed data {
  int N_latent = 750;
  
  int L1 = 5;
  int L2 = 25;
  int L3 = 125;
  
  vector[L1] factor1_level_probs = rep_vector(1.0 / L1, L1);
  vector[L2] factor2_level_probs = rep_vector(1.0 / L2, L2);
  vector[L3] factor3_level_probs = rep_vector(1.0 / L3, L3);
  
  real baseline = 5;
  real delta_factor1_tau = 2;
  
  real delta_factor1[L1] 
    = normal_rng(rep_vector(0, L1), delta_factor1_tau);

  real<lower=0> delta_factor2_tau[L1]
    = fabs(normal_rng(rep_vector(0, L1), 2.25));
    
  real<lower=0> delta_factor3_tau[L2]
      = fabs(normal_rng(rep_vector(0, L2), 1.5));
}

generated quantities {
  real sigma = 1;
    
  int<lower=1> N = N_latent;
  int<lower=1> N_factor1_levels = L1;
  int<lower=1> N_factor2_levels = L2;
  int<lower=1> N_factor3_levels = L3;
  
  int<lower=1, upper=N_factor3_levels> obs_to_factor3[N_latent];
  int<lower=1, upper=N_factor2_levels> factor3_to_factor2[L3];
  int<lower=1, upper=N_factor1_levels> factor2_to_factor1[L2];

  vector[N_latent] y;  

  { 
    vector[L2] delta_factor2;
    vector[L3] delta_factor3;
    
    for (l2 in 1:L2) {
      int l1 = categorical_rng(factor1_level_probs);
      factor2_to_factor1[l2] = l1;
      delta_factor2[l2] = normal_rng(0, delta_factor2_tau[l1]);
    }
    
    for (l3 in 1:L3) {
      int l2 = categorical_rng(factor2_level_probs);
      factor3_to_factor2[l3] = l2;
      delta_factor3[l3] = normal_rng(0, delta_factor3_tau[l2]);
    }

    for (n in 1:N) {
      int l3 = categorical_rng(factor3_level_probs);
      int l2 = factor3_to_factor2[l3];
      int l1 = factor2_to_factor1[l2];
      obs_to_factor3[n] = l3;

      y[n] = normal_rng(  baseline 
                        + delta_factor1[l1] 
                        + delta_factor2[l2] 
                        + delta_factor3[l3], sigma);
    }   
  }
}
