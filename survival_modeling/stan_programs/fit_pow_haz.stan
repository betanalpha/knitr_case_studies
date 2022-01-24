functions {
  // Hazard function
  real log_hazard(real t, real alpha, real sigma) {
    return log(alpha) - log(sigma) + (alpha - 1) * log(t / sigma);
  }
  
  // Survival function
  real log_survival(real t, real alpha, real sigma) {
    return -pow(t / sigma, alpha);
  }
  
  // Inverse survival function
  real inv_survival(real u, real alpha, real sigma) {
    return sigma * pow(-log(u), 1 / alpha);
  }
}

data {
  int<lower=1> N;   // Number of observations
  real obs_times[N]; // Observation times
  
  int<lower=1> N_display;                 // Number of display times
  real<lower=0> display_times[N_display]; // Display times
}

parameters {
  real<lower=0> alpha;
  real<lower=0> sigma;
}

model {
  alpha ~ normal(0, 5);
  sigma ~ normal(0, 5);
  
  for (n in 1:N)
    target +=   log_survival(obs_times[n], alpha, sigma)
              + log_hazard(obs_times[n], alpha, sigma);
}

generated quantities {
  real display_survival[N_display];
  real<lower=0> pred_obs_times[N];
  
  for (n in 1:N_display)
    display_survival[n] = exp(log_survival(display_times[n], alpha, sigma));
    
  for (n in 1:N)
    pred_obs_times[n] = inv_survival(uniform_rng(0, 1), alpha, sigma);
}
