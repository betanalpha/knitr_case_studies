data {
  int<lower=1> N_predict;
  real f[N_predict];
  real<lower=0> sigma;
}

parameters {}
model {}

generated quantities {
  real[N_predict] y_predict normal_rng(f, sigma);
}
