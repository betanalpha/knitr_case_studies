functions {
  real t_fall(real g, real h, real v0) {
    return (v0 + sqrt(square(v0) + 2 * g * h)) / g;
  }
}

data {
  // Detector calibration data
  int N_calibrations;
  real obs_t_calib[N_calibrations];    // seconds
  
  // Falling data
  int N_heights;
  real<lower=0> heights[N_heights];    // meters
  int N_drops;
  real obs_t_fall[N_drops, N_heights]; // seconds
}

parameters {
  real<lower=0> sigma_t; // seconds
  real v0[N_drops];      // meters per second
  real<lower=0> sigma_v; // meters per second
  real g;                // meters per second squared
}

model {
  // Prior model
  sigma_t ~ normal(0, 0.5);
  sigma_v ~ normal(0, 1);
  v0 ~ normal(0, sigma_v);
  g ~ inv_gamma(1.49, 2.83);

  // Observational model
  obs_t_calib ~ normal(rep_vector(0, N_calibrations), sigma_t);
  
  for (n in 1:N_drops) {
    for (h in 1:N_heights) {
      real mu_t_fall = t_fall(g, heights[h], v0[n]);
      obs_t_fall[n, h] ~ normal(mu_t_fall, sigma_t);
    }
  }
}

// Simulate a full observation from each value of the parameters
generated quantities {
  real obs_t_fall_post_pred[N_drops, N_heights];

  for (n in 1:N_drops) {
    for (h in 1:N_heights) {
      real mu_t_fall = t_fall(g, heights[h], v0[n]);
      obs_t_fall_post_pred[n, h] = normal_rng(mu_t_fall, sigma_t);
    }
  }
}
