data {
  int<lower=0> N;      // Number of observations
  vector[N] y;         // Observations
  real<lower=0> sigma; // Measurement variability
  
  // Number of individual contexts in hierarchy
  int<lower=0> K;
  // Individual context from which each observation is generated
  int<lower=1, upper=K> indiv_idx[N]; 
}

parameters {
  real mu;            // Population location
  real<lower=0> tau;  // Population scale
  vector[K] eta;      // Non-centered individual parameters
}

transformed parameters {
  // Recentered individual parameters
  vector[K] theta = mu + tau * eta;
}

model {
  mu ~ normal(0, 5);                   // Prior model
  tau ~ normal(0, 5);                  // Prior model
  eta ~ normal(0, 1);                  // Non-centered hierarchical model
  y ~ normal(theta[indiv_idx], sigma); // Observational model
}

// Simulate a full observation from the current value of the parameters
generated quantities {
  real y_post_pred[N] = normal_rng(theta[indiv_idx], sigma);
}
