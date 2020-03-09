data {
  int<lower=1> N; // Number of ys, supplied by interface
}

parameters {
  real y[N];  // Probabilistic variables y
  real theta; // Probabilistic variable theta
}

model {
  // Add conditional probability density for the ys
  // given theta to the joint log probabilty density
  // using the vectorized log probabilty density
  target += normal_lpdf(y | theta, 1);
  
  // Add marginal probabiltiy density for theta
  // to the joint log probability density
  target += normal_lpdf(theta | 0, 1);
}
