data {
  int<lower=0> N_obs; // Number of sales observations
  int y_obs[N_obs];  // Daily units sold
  
  int<lower=0> N_obs_days;                           // Number of days observed
  int<lower=1, upper=N_obs_days> obs_day_idx[N_obs]; // Day of each observation
  
  int<lower=0> N_pred_days;                        // Number of future days to predict
  int<lower=N_obs_days> pred_day_idx[N_pred_days]; // Day of each prediction

  vector<lower=0, upper=1>[N_obs_days] w; // Partial centering parameters
}

parameters {
  //vector[K_cp] delta_days_cp;   // Daily excess
  //vector[K_ncp] delta_days_ncp; // Daily excess
  vector[N_obs_days] delta_days_tilde; // Latent parameter for daily excess
  
  real<lower=0> inner_tau_day;              // Inner weekly log variation scale
  real<lower=inner_tau_day> outer_tau_day;  // Outer weekly log variation scale
  real<lower=0, upper=1> gamma; // Inner core probability
  
  real alpha;                   // Baseline sales per day
  real<lower=0> delta_season;   // Amplitude of log seasonal variation
  real delta_phi;               // Phase of seasonal variation from 10.5 days
  real<lower=0> milli_inv_psi;  // Observational dispersion times 1000
}

transformed parameters {
  real phi = 10.5 + delta_phi;
  
  vector[N_obs_days] delta_days;
  for (d in 1:N_obs_days)
    delta_days[d] = pow(inner_tau_day, 1 - w[d]) * delta_days_tilde[d];
}

model {
  for (d in 1:N_obs_days) {
    target += log_mix(gamma,
                      normal_lpdf(delta_days_tilde[d] | 0, pow(inner_tau_day, w[d])), 
                      normal_lpdf(delta_days_tilde[d] | 0, pow(inner_tau_day, w[d] - 1) * outer_tau_day));
  }

  inner_tau_day ~ normal(0, log(1.1));
  outer_tau_day ~ normal(0, log(1.5));
  gamma ~ beta(5, 0.5); 

  alpha ~ normal(log(200), log(50));
  delta_season ~ normal(0, log(1.5));
  delta_phi ~ normal(0, 1);
  milli_inv_psi ~ normal(0, 10);
  
  // Observational model
  y_obs 
    ~ neg_binomial_2_log(  alpha
                         + delta_season * sin(2 * pi() * to_vector(obs_day_idx) / 360 + phi)
                         + delta_days[obs_day_idx], 1000 / milli_inv_psi);
}

generated quantities {
  real log_mu[N_obs_days];
  int y_post_retro[N_obs_days];
  int y_post_pre[N_pred_days];

  for (d in 1:N_obs_days) {
    real lm = alpha + delta_season * sin(2 * pi() * d / 360 + phi) + delta_days[d];
    
    // Truncate extreme samples to avoid numerical overflow 
    // in neg_binomial_2_log pseudo random number generator
    if (lm > +20) lm = +20;
    if (lm < -20) lm = -20;
    
    log_mu[d] = lm;
    y_post_retro[d] = neg_binomial_2_log_rng(lm, 1000 / milli_inv_psi);
  }
  
  for (d in 1:N_pred_days) {
    real delta_day = 0;
    real lm = 0;
    
    if (bernoulli_rng(gamma))
      delta_day = normal_rng(0, inner_tau_day);
    else 
      delta_day = normal_rng(0, outer_tau_day);
    
    lm =   alpha 
         + delta_season * sin(2 * pi() * pred_day_idx[d] / 360 + phi)
         + delta_day;
              
    // Truncate extreme samples to avoid numerical overflow 
    // in neg_binomial_2_log pseudo random number generator
    if (lm > +20) lm = +20;
    if (lm < -20) lm = -20;
    
    y_post_pre[d] = neg_binomial_2_log_rng(lm, 1000 / milli_inv_psi);
  }
}
