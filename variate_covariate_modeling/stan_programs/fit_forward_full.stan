data {
  int<lower=1> N;
  vector[N] x1; // Complete observation
  vector[N] y1; // Complete observation
  vector[N] x2; // Incomplete observation
}

parameters { 
  real mu1;
  real<lower=0> tau1;
  
  real mu2;
  real<lower=0> tau2;
  
  real<lower=0> sigma;
}

model {
  mu1 ~ normal(0, 3);
  tau1 ~ normal(0, 1);
  
  mu2 ~ normal(0, 3);
  tau2 ~ normal(0, 1);
  
  sigma ~ normal(0, 1);
  
  x1 ~ normal(mu1, tau1);
  y1 ~ normal(x1, sigma);
  
  x2 ~ normal(mu2, tau2);
}

generated quantities {
  // Predict completion of second observation
  real y2_pred[N] = normal_rng(x2, sigma);
}
