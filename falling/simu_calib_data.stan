data {
  int N_calibrations;
}

transformed data {
  real sigma_t = 0.02; // seconds
}

generated quantities {
  real obs_t[N_calibrations];
  for (n in 1:N_calibrations)
      obs_t[n] = normal_rng(0, sigma_t);
}
