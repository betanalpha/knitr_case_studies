data {
 int<lower = 0> N;
 vector[N] y;
}

transformed data {
  real<lower=0> sigma[2] = {0.75, 1.25};
}

parameters {
  vector[2] mu;
  real<lower=0, upper=1> theta;
}

model {
  mu[1] ~ normal(-2, 2);
  mu[2] ~ normal(2, 2);
  theta ~ beta(3, 3);
  for (n in 1:N)
    target += log_mix(theta,
                      normal_lpdf(y[n] | mu[1], sigma[1]),
                      normal_lpdf(y[n] | mu[2], sigma[2]));
}
