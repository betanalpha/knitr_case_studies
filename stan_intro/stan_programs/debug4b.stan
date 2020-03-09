transformed data {
  int N = 10;
  real mu[N];
  for (n in 1:N)
    mu[N] = n;
  print("mu = ", mu);
}

parameters {
  real x[N];
}

model {
  x ~ normal(mu, 1);
}
