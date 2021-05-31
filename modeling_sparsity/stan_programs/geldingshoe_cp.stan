data {
  int<lower=0> K;                       // Number of groups
  int<lower=0> N;                       // Number of observations
  vector[N] y;                          // Observations
  int<lower=1, upper=K> context_idx[N]; // Context assignments
  real<lower=0> tau;                    // Base Population Scale
  real<lower=0> sigma;                  // Measurement variability
}

parameters {
  vector[K] theta;
  vector<lower=0>[K] lambda;
}

model {
  // Horseshoe prior model
  theta ~ normal(0, lambda * tau);
  lambda ~ cauchy(0, 1);

  // Observational model
  y ~ normal(theta[context_idx], sigma);
}
