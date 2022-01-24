functions {
  // Inverse survival function
  real inv_survival(real u, real alpha, real sigma) {
    return sigma * pow(-log(u), 1 / alpha);
  }
}

data {
  int<lower=1> N;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

generated quantities {
  real<lower=0> obs_times[N];
  
  for (n in 1:N) {
    real u = uniform_rng(0, 1);
    obs_times[n] = inv_survival(u, alpha, sigma);
  }
}
