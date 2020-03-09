transformed data {
  int N = 10;
  real mu[N];
  for (n in 1:N)
    mu[N] = n;
}

parameters {
  real x[N];
}

model {
  x ~ normal(mu, 1);
}
