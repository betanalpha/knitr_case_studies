data {
  int<lower=1> N;
  real t[N + 1];
  
  int<lower=1> N_obs;
  int<lower=1, upper=N + 1> obs_idx[N_obs];
  int<lower=1, upper=N + 1> unobs_idx[N - N_obs];
  real y_obs[N_obs];
  real<lower=0> sigma;
}

parameters {
  real mu;
  real<lower=0> gamma;
  real<lower=0> tau;
  
  real x0;
  real x_cp[N_obs];
  real x_ncp[N - N_obs];
}

transformed parameters {
  real x[N + 1]; 
  
  x[1] = x0;
  x[obs_idx] = x_cp;
  
  for (n in 1:(N - N_obs)) {
    int m = unobs_idx[n];  
    real epsilon = t[m] - t[m - 1];
    x[m] =  x[m - 1] + gamma * (mu - x[m - 1]) * epsilon 
          + tau * sqrt(epsilon) * x_ncp[n];
  }
}

model {
  mu ~ normal(0, 10);
  gamma ~ normal(0, 1);
  tau ~ normal(0, 1);
  
  x0 ~ normal(0, 3);
  
  for (n in 1:N_obs) {
    int m = obs_idx[n];
    real epsilon = t[m] - t[m - 1];
    x_cp[n] ~ normal(x[m - 1] + gamma * (mu - x[m - 1]) * epsilon, 
                     tau * sqrt(epsilon) );
  }
  
  x_ncp ~ normal(0, 1);

  y_obs ~ normal(x[obs_idx], sigma);
}
