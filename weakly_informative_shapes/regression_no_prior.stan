data {
  int<lower=1> N;
  vector[N] x; // Rainfall in cm
  vector[N] y; // Income in k$
}

parameters {
  real alpha;          // k$
  real beta;           // k$ / cm
  real<lower=0> sigma; // k$
}

model {
  y ~ normal(beta * x + alpha, sigma);
}
