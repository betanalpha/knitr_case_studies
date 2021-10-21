data {
  real t0;
  
  int<lower=1> N_obs;
  real t_obs[N_obs];
  real y_obs[N_obs];
  real<lower=0> sigma;
}

parameters {
  real mu;
  real<lower=0> gamma;
  real<lower=0> tau;

  real x0;
  real x[N_obs];
}

model {
  real decay;
  
  mu ~ normal(0, 10);
  gamma ~ normal(0, 1);
  tau ~ normal(0, 1);
  
  x0 ~ normal(0, 3);
  
  decay = exp(-gamma * (t_obs[1] - t0));
  x[1] ~ normal(mu + (x0 - mu) * decay, 0.7071 * tau * sqrt( (1 - decay) / gamma ) );
  
  for (n in 2:N_obs) {
    decay = exp(-gamma * (t_obs[n] - t_obs[n - 1]));
    x[n] ~ normal(mu + (x[n - 1] - mu) * decay, 0.7071 * tau * sqrt( (1 - decay) / gamma ) );
  }
  
  y_obs ~ normal(x, sigma);
}
