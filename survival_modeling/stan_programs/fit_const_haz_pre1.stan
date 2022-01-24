functions {
  // Hazard function
  real log_const_haz_hazard(real t, real gamma) {
    return log(gamma);
  }
  
  // Survival function
  real log_const_haz_survival(real t, real gamma) {
    return -gamma * t;
  }
  
  // Inverse survival function
  real inv_const_haz_survival(real u, real gamma) {
    return - log(u) / gamma;
  }
}

data {
  int<lower=1> N;             // Number of observations
  real<lower=0> obs_times[N]; // Observation times
}

parameters {
  real<lower=0> gamma1;
  real<lower=0> gamma2;
  
  real<lower=0, upper=1> t1_tilde[N];
}

transformed parameters {
  // Implement varying constraints
  real t1[N];
  for (n in 1:N) t1[n] = t1_tilde[n] * obs_times[n];
}

model {
  gamma1 ~ normal(0, 5);
  gamma2 ~ normal(0, 5);
  
  for (n in 1:N) {
    // Unobserved precursor event
    target +=   log_const_haz_survival(t1[n], gamma1)
              + log_const_haz_hazard(t1[n], gamma1);
              
    // Observed event
    target +=   log_const_haz_survival(obs_times[n] - t1[n], gamma2)
              + log_const_haz_hazard(obs_times[n] - t1[n], gamma2);
  }
}

generated quantities {
  real<lower=0> pred_obs_times[N];
  
  for (n in 1:N) {
    real t = inv_const_haz_survival(uniform_rng(0, 1), gamma1);
    pred_obs_times[n] = inv_const_haz_survival(uniform_rng(0, 1), gamma2) + t;
  }
}
