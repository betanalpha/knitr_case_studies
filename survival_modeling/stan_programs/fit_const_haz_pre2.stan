functions {
  // Inverse survival function
  real inv_const_haz_survival(real u, real gamma) {
    return - log(u) / gamma;
  }
  
  // Observational model for subsequent events given unobserved precursor events
  real sub_lpdf(real t, real gamma1, real gamma2) {
    real diff = gamma2 - gamma1;
    if (fabs(diff) < 0.0001) {
      return log(gamma1) + log(gamma2) - gamma2 * t + log(t);
    } else if(diff > 0) {
      return   log(gamma1) + log(gamma2) 
             + log_diff_exp(-gamma1 * t, -gamma2 * t) 
             - log(diff);
    } else {
      return   log(gamma1) + log(gamma2) 
             + log_diff_exp(-gamma2 * t, -gamma1 * t) 
             - log(-diff);
    }
  }
  
  // Survival probability for subsequent events
  real sub_survival(real t, real gamma1, real gamma2) {
     return   ( gamma2 * exp(-gamma1 * t) - gamma1 * exp(-gamma2 * t) ) 
            / ( gamma2 - gamma1 );
  }
}

data {
  int<lower=1> N;             // Number of observations
  real<lower=0> obs_times[N]; // Observation times
  
  int<lower=1> N_display;                 // Number of display times
  real<lower=0> display_times[N_display]; // Display times
}

parameters {
  real<lower=0> gamma1;
  real<lower=0> gamma2;
}

model {
  gamma1 ~ normal(0, 5);
  gamma2 ~ normal(0, 5);
  
  for (n in 1:N)
    target += sub_lpdf(obs_times[n] | gamma1, gamma2);
}

generated quantities {
  real display_survival[N_display];
  real<lower=0> pred_obs_times[N];
  
  for (n in 1:N_display)
    display_survival[n] = sub_survival(display_times[n], gamma1, gamma2);
    
  for (n in 1:N) {
    real t = inv_const_haz_survival(uniform_rng(0, 1), gamma1);
    pred_obs_times[n] = inv_const_haz_survival(uniform_rng(0, 1), gamma2) + t;
  }
}
