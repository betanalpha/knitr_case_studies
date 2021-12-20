data {
  int<lower=1> N;
  vector[N] x1; // Complete observation
  vector[N] y1; // Complete observation
  vector[N] x2; // Incomplete observation
}

parameters { 
  real<lower=0> sigma;
}

model {
  sigma ~ normal(0, 1);
  
  y1 ~ normal(x1, sigma);
}

generated quantities {
  // Predict completion of second observation using
  // conditional variate model from first observation
  real y2_pred[N] = normal_rng(x2, sigma);
}
