data {
  int N_calibrations;
  real obs_t_calib[N_calibrations]; // seconds
}

parameters {
  real<lower=0> sigma_t; // seconds
}

model {
  // Prior model
  sigma_t ~ normal(0, 0.5);

  // Observational model
  obs_t_calib ~ normal(rep_vector(0, N_calibrations), sigma_t);
}
