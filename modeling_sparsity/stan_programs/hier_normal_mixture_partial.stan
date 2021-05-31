data {
  int<lower=0> K;                       // Number of groups
  int<lower=0> N;                       // Number of observations
  vector[N] y;                          // Observations
  int<lower=1, upper=K> context_idx[N]; // Context assignments
  vector<lower=0, upper=1>[K] w;        // Partial centering parameters

  real<lower=0> sigma;                  // Measurement variability
}

parameters {
  vector[K] theta_tilde;
  
  real<lower=0, upper=1> gamma;
  real<lower=0> tau1;
  real<lower=tau1> tau2;
}

transformed parameters {
  vector[K] theta;
  for (k in 1:K)
    theta[k] = pow(tau1, 1 - w[k]) * theta_tilde[k];
}

model {
  for (k in 1:K) {
    target += log_mix(gamma,
                      normal_lpdf(theta_tilde[k] | 0, pow(tau1, w[k])), 
                      normal_lpdf(theta_tilde[k] | 0, pow(tau1, w[k] - 1) * tau2));
  }
  
  gamma ~ beta(3, 1);
  tau1 ~ normal(0, 0.1);
  tau2 ~ normal(0, 10);

  // Observational model
  y ~ normal(theta[context_idx], sigma);
}
