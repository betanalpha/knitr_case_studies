// Values provided externally
data {
  real mu;
  real<lower=0> sigma;
}

generated quantities {
  // Generate an exact sample from a distribution
  // specified by a normal probabilty density function
  real x = normal_rng(mu, sigma);
}
