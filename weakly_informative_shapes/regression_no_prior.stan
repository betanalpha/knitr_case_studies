data {
  real x; // Rainfall in cm
  real y; // Income in k$
}

parameters {
  real alpha;          // k$
  real beta;           // k$ / cm
  real<lower=0> sigma; // k$
}

model {
  y ~ normal(beta * x + alpha, sigma);
}
