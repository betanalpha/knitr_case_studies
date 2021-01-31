transformed data {
  int N_latent = 300;
  
  int N1 = 5;
  int N2 = 8;
  int N3 = 4;
  
  real baseline = 5;
  
  real<lower=0> delta1_tau = 4;
  vector[N1] delta1 = to_vector(normal_rng(rep_vector(0, N1), delta1_tau));

  real<lower=0> delta2_tau = 7;
  vector[N2] delta2 = to_vector(normal_rng(rep_vector(0, N2), delta2_tau));

  real<lower=0> delta3_tau = 5;
  vector[N3] delta3 = to_vector(normal_rng(rep_vector(0, N3), delta3_tau));

  // First-order interactions
  int<lower=0> N_factor12_levels = N1 * N2;

  real<lower=0> delta12_tau = 3;
  vector[N1 * N2] delta12 = to_vector(normal_rng(rep_vector(0, N1 * N2), delta12_tau));

  int<lower=0> N_factor13_levels = N1 * N3;

  real<lower=0> delta13_tau = 0.25;
  vector[N1 * N3] delta13 = to_vector(normal_rng(rep_vector(0, N1 * N3), delta13_tau));

  int<lower=0> N_factor23_levels = N2 * N3;

  real<lower=0> delta23_tau = 2;
  vector[N2 * N3] delta23 = to_vector(normal_rng(rep_vector(0, N2 * N3), delta23_tau));

  // Second-order interactions
  int<lower=0> N_factor123_levels = N1 * N2 * N3;
  
  real<lower=0> delta123_tau = 0.01;
  vector[N_factor123_levels] delta123 = to_vector(normal_rng(rep_vector(0, N1 * N2 * N3), delta123_tau));
  
  vector[N1] factor1_level_probs = rep_vector(1.0, N1);
  vector[N2] factor2_level_probs = rep_vector(1.0, N2);
  vector[N3] factor3_level_probs = rep_vector(1.0, N3);
  
  factor1_level_probs[1] = 6;
  factor1_level_probs[4] = 3;
  factor1_level_probs /= sum(factor1_level_probs);
  
  factor2_level_probs[1] = 2;
  factor2_level_probs[3] = 5;
  factor2_level_probs[6] = 4;
  factor2_level_probs /= sum(factor2_level_probs);
  
  factor3_level_probs[2] = 3;
  factor3_level_probs /= sum(factor3_level_probs);
}

generated quantities {
  int N = N_latent;
  
  int<lower=0> N_main_factors = 3;
  int<lower=0> N_main_factor_levels[3] = {N1, N2, N3};
  int<lower=1> main_factor_level_idx[3, N_latent];
    
  real y[N_latent];
  real sigma = 1;

  { 
    vector[N_latent] theta;
    
    
    int factor1_idx[N_latent];
    int factor2_idx[N_latent];
    int factor3_idx[N_latent];
    
    int factor12_idx[N_latent];
    int factor13_idx[N_latent];
    int factor23_idx[N_latent];
    int factor123_idx[N_latent];

    for (n in 1:N) {
      factor1_idx[n] = categorical_rng(factor1_level_probs);
      factor2_idx[n] = categorical_rng(factor2_level_probs);
      factor3_idx[n] = categorical_rng(factor3_level_probs);
      
      factor12_idx[n] = (factor1_idx[n] - 1) * N2 + factor2_idx[n];
      factor23_idx[n] = (factor2_idx[n] - 1) * N3 + factor3_idx[n];
      factor13_idx[n] = (factor1_idx[n] - 1) * N3 + factor3_idx[n];
     
      factor123_idx[n] =   (factor1_idx[n] - 1) * N2 * N3
                        + (factor2_idx[n] - 1) * N3 
                        +  factor3_idx[n];
    }
      
    main_factor_level_idx[1] = factor1_idx;
    main_factor_level_idx[2] = factor2_idx;
    main_factor_level_idx[3] = factor3_idx;

    theta =  baseline
           + delta1[factor1_idx]
           + delta2[factor2_idx]
           + delta3[factor3_idx]
           + delta12[factor12_idx]
           + delta13[factor13_idx]
           + delta23[factor23_idx]
           + delta123[factor123_idx];
                      
    y = normal_rng(theta, sigma);   
  }
}
