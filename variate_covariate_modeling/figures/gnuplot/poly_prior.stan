data {
  int<lower=1> N;
  vector[N] x;
  vector[N] y;

  int<lower=1> N_predict;
  real x_predict[N_predict];
  
  real<lower=0> sigma;
}

parameters {
  real c0;
  real c1;
  real c2;
  real c3;
}

model {
  c0 ~ normal(0, 10);
  c1 ~ normal(0, 2);
  c2 ~ normal(0, 0.01);
  c3 ~ normal(0, 0.01);
}

generated quantities {
  vector[N_predict] f;
  for (n in 1:N_predict) {
    f[n] =   c0
           + c1 * x_predict[n]
           + c2 * x_predict[n] * x_predict[n]
           + c3 * x_predict[n] * x_predict[n] * x_predict[n];
  }
}
