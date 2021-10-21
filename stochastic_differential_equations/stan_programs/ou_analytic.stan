data {
  int<lower=1> N;
  real t[N + 1];
  
  real x0;

  real mu;
  real<lower=0> gamma;
  real<lower=0> tau;
}

parameters {
  real x_free[N];
}

transformed parameters {
  real x[N + 1] = append_array(rep_array(x0, 1), x_free);
}

model {
  for (n in 1:N) {
    real decay = exp(-gamma * (t[n + 1] - t[n]));
    x[n + 1] ~ normal(mu + (x[n] - mu) * decay, 0.7071 * tau * sqrt( (1 - decay) / gamma ) );
  }
}
