transformed data {
  int N;
  real mu[N];
  for (n in 1:N)
    mu[n] = n;
}

parameters {
  real x[N];
}

model {
  x ~ normal(mu, 1);
}
