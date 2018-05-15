data {
  int N;
  int y[N];
}

parameters {
  real<lower=0> lambda;
}

model {
  lambda ~ normal(0, 6.44787);
  y ~ poisson(lambda);
}
