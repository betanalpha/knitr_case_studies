data {
 int<lower = 0> N;
 vector[N] y;
}

parameters {
  vector[2] mu;
  real<lower=0> sigma[2];
  real<lower=0, upper=1> theta;
}

model {
 sigma ~ normal(0, 2);
 mu[1] ~ normal(4, 0.5);
 mu[2] ~ normal(-4, 0.5);
 theta ~ beta(5, 5);
 for (n in 1:N)
   target += log_mix(theta,
                     normal_lpdf(y[n] | mu[1], sigma[1]),
                     normal_lpdf(y[n] | mu[2], sigma[2]));
}
