parameters {
  real a;          // Unknown truncation
  real<lower=a> x; // Constrained parameter
}

model {
  target +=   normal_lpdf(a | 0, 1);  // Truncation density
  target +=   normal_lpdf(x | 0, 1);  // Nominal target density
  target += - normal_lccdf(a | 0, 1); // Corrected normalization
}
