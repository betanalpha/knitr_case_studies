functions {
  // Inverse survival function
  real inv_survival(real u, real gamma) {
    return - log(u) / gamma;
  }
}

data {
  int<lower=1> N;
  real<lower=0> gamma1;
  real<lower=0> gamma2;
}

generated quantities {
  real<lower=0> obs_times[N];
  
  for (n in 1:N) {
    // Simulate precursor event
    real u1 = uniform_rng(0, 1);
    real t1 = inv_survival(u1, gamma1);
    // Simulate subsequent event conditional on precursor event
    real u2 = uniform_rng(0, 1);
    obs_times[n] = inv_survival(u2, gamma2) + t1;
  }
}
