data {
  int<lower=1> N;
  real f[N];
  real<lower=0> sigma;
}

parameters {}
model {}

generated quantities {
  vector[N] y_predict;
  for (n in 1:N)
    y_predict[n] = normal_rng(f[n], sigma);
}
