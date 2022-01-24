functions {
  // Survival function
  real survival(real t, real gamma) {
    return exp(-gamma * t);
  }
  
  // Inverse survival function
  real inv_survival(real u, real gamma) {
    return - log(u) / gamma;
  }
}

data {
  int<lower=1> N2;
  real<lower=0> gamma;
  
  real t1; // Left-censoring threshold
  real t2; // Interval-censoring left threshold
  real t3; // Interval-censoring right threshold/Right-censoring threshold
}

transformed data {
  real p1 = 1                   - survival(t1, gamma);
  real p2 = survival(t1, gamma) - survival(t2, gamma);
  real p3 = survival(t2, gamma) - survival(t3, gamma);
  real p4 = survival(t3, gamma) - 0;

  // Be careful about Stan's non-standard negative binomial parameterization
  int N_remaining = neg_binomial_rng(N2, p2 / (1 - p2));
  
  vector[3] cens_probs = [p1 / (1 - p2), p3 / (1 - p2), p4 / (1 - p2)]';
  int Ns[3] = multinomial_rng(cens_probs, N_remaining);
}

generated quantities {
  int<lower=0> N1 = Ns[1];
  real<lower=0> obs_times[N2] = rep_array(0.0, N2);
  int<lower=0> N3 = Ns[2];
  int<lower=0> N4 = Ns[3];
  
  for (n in 1:N2) {
    real u = uniform_rng(survival(t2, gamma), survival(t1, gamma));
    obs_times[n] = inv_survival(u, gamma);
  }
}
