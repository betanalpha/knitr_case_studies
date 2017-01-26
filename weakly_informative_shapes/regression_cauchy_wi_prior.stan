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
  alpha ~ cauchy(0, 1);
  beta ~ cauchy(0, 1);
  sigma ~ cauchy(0, 1);
  y ~ normal(beta * x + alpha, sigma);
}
