functions {
  // Hazard function
  real log_hazard(real t, real gamma) {
    return log(gamma);
  }
  
  // Survival function
  real log_survival(real t, real gamma) {
    return -gamma * t;
  }
  
  // Inverse survival function
  real inv_survival(real u, real gamma) {
    return - log(u) / gamma;
  }
}

data {
  int<lower=1> N;             // Number of observations
  real<lower=0> obs_times[N]; // Observation times
  
  int<lower=1> N_display;                 // Number of display times
  real<lower=0> display_times[N_display]; // Display times
}

parameters {
  real<lower=0> gamma;
}

model {
  gamma ~ normal(0, 5);
  
  for (n in 1:N)
    target +=   log_survival(obs_times[n], gamma)
              + log_hazard(obs_times[n], gamma);
}

generated quantities {
  real display_survival[N_display];
  real<lower=0> pred_obs_times[N];
  
  for (n in 1:N_display)
    display_survival[n] = exp(log_survival(display_times[n], gamma));
    
  for (n in 1:N)
    pred_obs_times[n] = inv_survival(uniform_rng(0, 1), gamma);
}
