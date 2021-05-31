data {
  int<lower=0> K;                       // Number of groups
  int<lower=0> N;                       // Number of observations
  vector[N] y;                          // Observations
  int<lower=1, upper=K> context_idx[N]; // Context assignments
  
  vector<lower=0, upper=1>[K] w;        // Partial centering parameters
  
  real<lower=0> tau;                    // Base Population Scale
  real<lower=0> sigma;                  // Measurement variability
}

parameters {
  vector[K] theta_tilde;
  vector<lower=0>[K] lambda;
}

transformed parameters {
  vector[K] theta;
  for (k in 1:K)
    theta[k] = pow(lambda[k] * tau, 1 - w[k]) * theta_tilde[k];
}

model {
  for (k in 1:K)
    theta_tilde[k] ~ normal(0, pow(lambda[k] * tau, w[k]));
  lambda ~ cauchy(0, 1);

  // Observational model
  y ~ normal(theta[context_idx], sigma);
}
