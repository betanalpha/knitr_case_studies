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
  int<lower=1> N1; // Number of observations in left-censored interval
  int<lower=1> N2; // Number of observations in uncensored interval
  int<lower=1> N3; // Number of observations in censored interval
  int<lower=1> N4; // Number of observations in right-censored interval
      
  real<lower=0> obs_times[N2]; // Observation times
  
  real t_min;
  real<lower=t_min> t1; // Left-censoring threshold
  real<lower=t1>    t2; // Interval-censoring left threshold
  real<lower=t2>    t3; // Interval-censoring right threshold/Right-censoring threshold
  real<lower=t3> t_max;
  
  int<lower=1> N_display;                 // Number of display times
  real<lower=0> display_times[N_display]; // Display times
}

transformed data {
  int N = N1 + N2 + N3 + N4;
}

parameters {
  real<lower=0> gamma;
}

model {
  real log_S1 = log_survival(t1, gamma);
  real log_S2 = log_survival(t2, gamma);
  real log_S3 = log_survival(t3, gamma);
  
  real log_p1 = log_diff_exp(0, log_S1);
  real log_p3 = log_diff_exp(log_S2, log_S3);
  real log_p4 = log_S3;
  
  gamma ~ normal(0, 5);
  
  target += N1 * log_p1;
  
  for (n in 1:N2)
    target +=   log_survival(obs_times[n], gamma)
              + log_hazard(obs_times[n], gamma);
              
  target += N3 * log_p3;
  target += N4 * log_p4;
}

generated quantities {
  real display_survival[N_display];
  real<lower=0> pred_obs_times[N];
  
  for (n in 1:N_display)
    display_survival[n] = exp(log_survival(display_times[n], gamma));
    
  for (n in 1:N) {
    real u = uniform_rng(0, 1);
    real time = inv_survival(u, gamma);
      
    if (time < t1) {
      pred_obs_times[n] = 0.5 * (t_min + t1);
    } else if (time < t2) {
      pred_obs_times[n] = time;
    } else if (time < t3) {
      pred_obs_times[n] = 0.5 * (t2 + t3);
    } else {
      pred_obs_times[n] = 0.5 * (t3 + t_max);
    }
  }  
}
