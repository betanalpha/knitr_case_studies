functions {
  // Inverse survival function
  real inv_survival(real u, real gamma) {
    return - log(u) / gamma;
  }
}

data {
  int<lower=1> N;
  real<lower=0> gamma;
}

generated quantities {
  real<lower=0> obs_times[N];
  
  for (n in 1:N) {
    real u = uniform_rng(0, 1);
    obs_times[n] = inv_survival(u, gamma);
  }
}
