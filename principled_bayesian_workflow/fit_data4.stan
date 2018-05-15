data {
  int N;
  int y[N];
}

transformed data {
  int U = 14;
}

parameters {
  real<lower=0, upper=1> theta;
  real<lower=0> lambda;
}

model {
  theta ~ beta(2.8663, 2.8663);
  lambda ~ inv_gamma(3.48681, 9.21604);

  for (n in 1:N) {
    real lpdf = poisson_lpmf(y[n] | lambda) - poisson_lcdf(U | lambda);
    if (y[n] == 0)
      target += log_mix(theta, 0, lpdf);
    else
      target += log1m(theta) + lpdf;
  }
}
