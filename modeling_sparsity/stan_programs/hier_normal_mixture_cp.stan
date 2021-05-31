data {
  int<lower=0> K;                       // Number of groups
  int<lower=0> N;                       // Number of observations
  vector[N] y;                          // Observations
  int<lower=1, upper=K> context_idx[N]; // Context assignments
  real<lower=0> sigma;                  // Measurement variability
}

parameters {
  vector[K] theta;
  
  real<lower=0, upper=1> gamma;
  real<lower=0> tau1;
  real<lower=tau1> tau2;
}

model {
  for (k in 1:K)
    target += log_mix(gamma,
                      normal_lpdf(theta[k] | 0, tau1), 
                      normal_lpdf(theta[k] | 0, tau2));
  
  gamma ~ beta(3, 1);
  tau1 ~ normal(0, 0.1);
  tau2 ~ normal(0, 10);

  // Observational model
  y ~ normal(theta[context_idx], sigma);
}
