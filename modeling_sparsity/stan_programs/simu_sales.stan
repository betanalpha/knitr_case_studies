data {
  int<lower=0> N_obs; // Number of sales observations
  
  int<lower=0> N_obs_days;                           // Number of days observed
  int<lower=1, upper=N_obs_days> obs_day_idx[N_obs]; // Day of each observation
  
  int<lower=0> N_pred_days;                        // Number of future days to predict
  int<lower=N_obs_days> pred_day_idx[N_pred_days]; // Day of each prediction
}

transformed data {
  real alpha = log(226);                // Baseline sales per day
  real delta_season = log(1.3);         // Amplitude of seasonal variation
  real phi = 10;                        // Phase of seasonal variation
  real<lower=0> tau_day1 = log(1.005);  // Weekly variation inner scale
  real<lower=0> tau_day2 = log(1.25);   // Weekly variation outer scale
  real<lower=0, upper=1> gamma = 0.9;   // Sparsity proportion
  real<lower=0> inv_psi = 0.001;        // Observational dispersion
}

generated quantities {
  int y_obs[N_obs];
  int y_pred[N_pred_days];
  vector[N_obs_days] delta_days;

  for (d in 1:N_obs_days) {
    if (bernoulli_rng(gamma))
      delta_days[d] = normal_rng(0, tau_day1);
    else 
      delta_days[d] = normal_rng(0, tau_day2);
  }

  y_obs 
    = neg_binomial_2_log_rng(  alpha
                             + delta_season * sin(2 * pi() * to_vector(obs_day_idx) / 360 + phi)
                             + delta_days[obs_day_idx], 1 / inv_psi);

  for  (d in 1:N_pred_days) {
    real delta_day = 0;
    if (bernoulli_rng(gamma))
      delta_day = normal_rng(0, tau_day1);
    else 
      delta_day = normal_rng(0, tau_day2);

    y_pred[d] 
      = neg_binomial_2_log_rng(  alpha
                               + delta_season * sin(2 * pi() * pred_day_idx[d] / 360 + phi)
                               + delta_day, 1 / inv_psi);
  }
}
