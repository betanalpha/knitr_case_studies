data {
  int<lower=1> N; // Number of ys, supplied by interface
  real y[N]; // Observed values of the ys, supplied by the interface
}

parameters {
  real theta; // Probabilistic variable theta
}

model {
  // Add conditional probability density for the ys
  // given theta to the joint log probabilty density
  // using equivalent sampling statement
  y ~ normal(theta, 1);
  
  // Add marginal probabiltiy density for theta
  // to the joint log probability density using
  // equivalent sampling statement
  theta ~ normal(0, 1);
}

generated quantities {
  // Save the exponential of theta for analysis
  real exp_theta = exp(theta);
}
