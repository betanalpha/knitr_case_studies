transformed data {
  real a = 0.25; // Known truncation
}

parameters {
  real<lower=a> x; // Constrained parameter
}

model {
  target += normal_lpdf(x | 0, 1);  // Nominal target density
}
