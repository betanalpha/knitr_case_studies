data {
  int<lower=0> N;      // Number of observations
  vector[N] y;         // Observations
  real<lower=0> sigma; // Measurement variability
  
  // Factor information
  int<lower=0> N_main_factors;
  int<lower=0> N_main_factor_levels[N_main_factors];
  int<lower=1> main_factor_level_idx[N_main_factors, N];
  
  // First-order factor interaction contributions logistics
  int<lower=0> N_inter1_factors;
  int<lower=0> N_inter1_factor_levels[N_inter1_factors];
  int<lower=1> inter1_factor_level_idx[N_inter1_factors, N];

  // Second-order factor interaction contributions logistics
  int<lower=0> N_inter2_factors;
  int<lower=0> N_inter2_factor_levels[N_inter2_factors];
  int<lower=1> inter2_factor_level_idx[N_inter2_factors, N];
}

parameters {  
  // Zeroth-order contribution (Baseline)
  real baseline;
  
  // First-order contributions (Main Effects)
  vector[N_main_factor_levels[1]] delta1;
  real<lower=0> delta1_tau;
  
  vector[N_main_factor_levels[2]] delta2;
  real<lower=0> delta2_tau;
  
  vector[N_main_factor_levels[3]] delta3;
  real<lower=0> delta3_tau;
  
  // Second-order contributions (First-order interactions)
  vector[N_inter1_factor_levels[1]] delta12;
  real<lower=0> delta12_tau;
  
  vector[N_inter1_factor_levels[2]] delta13;
  real<lower=0> delta13_tau;
  
  vector[N_inter1_factor_levels[3]] delta23;
  real<lower=0> delta23_tau;
  
  // Third-order contributions (Second-order interactions)
  vector[N_inter2_factor_levels[1]] delta123;
  real<lower=0> delta123_tau;
}

model {
  vector[N] theta =  baseline
                   + delta1[main_factor_level_idx[1]]
                   + delta2[main_factor_level_idx[2]]
                   + delta3[main_factor_level_idx[3]]
                   + delta12[inter1_factor_level_idx[1]]
                   + delta13[inter1_factor_level_idx[2]]
                   + delta23[inter1_factor_level_idx[3]]
                   + delta123[inter2_factor_level_idx[1]];
    
  baseline ~ normal(0, 8);
  
  delta1 ~ normal(0, delta1_tau);
  delta1_tau ~ normal(0, 4);
  
  delta2 ~ normal(0, delta2_tau);
  delta2_tau ~ normal(0, 4);
  
  delta3 ~ normal(0, delta3_tau);
  delta3_tau ~ normal(0, 4);
  
  delta12 ~ normal(0, delta12_tau);
  delta12_tau ~ normal(0, 2);
  
  delta23 ~ normal(0, delta23_tau);
  delta23_tau ~ normal(0, 2);
  
  delta13 ~ normal(0, delta13_tau);
  delta13_tau ~ normal(0, 2);
  
  delta123 ~ normal(0, delta123_tau);
  delta123_tau ~ normal(0, 1);

  y ~ normal(theta, sigma);
}
