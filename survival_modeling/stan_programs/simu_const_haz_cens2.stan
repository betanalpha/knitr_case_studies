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
  int<lower=1> N;
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
  
  vector[4] probs = [p1, p2, p3, p4]';
  int Ns[4] = multinomial_rng(probs, N);
}

generated quantities {
  int<lower=0> N1 = Ns[1];
  int<lower=0> N2 = Ns[2];
  int<lower=0> N3 = Ns[3];
  int<lower=0> N4 = Ns[4];
  
  real<lower=0> obs_times[Ns[2]];
  
  for (n in 1:N2) {
    real u = uniform_rng(survival(t2, gamma), survival(t1, gamma));
    obs_times[n] = inv_survival(u, gamma);
  }
}
