parameters {
  real<lower=0> theta; // Constrained declaration for theta
}

model {
  // Gamma density valid only for positive values of theta
  theta ~ gamma(1.25, 1.25);
}
