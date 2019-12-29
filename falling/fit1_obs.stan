functions {
  real t_fall(real g, real h) {
    return sqrt(2 * h / g);
  }
}

data {
  int N_heights;
  real<lower=0> heights[N_heights];    // meters
  int N_drops;
  real obs_t_fall[N_drops, N_heights]; // seconds
}

parameters {
  real<lower=0> sigma_t; // seconds
  real g;                // meters per second squared
}

model {
  // Prior model
  sigma_t ~ ???;
  g ~ ???;

  // Observational model
  for (n in 1:N_drops) {
    for (h in 1:N_heights) {
      real mu_t_fall = t_fall_naive(g, heights[h]);
      obs_t_fall[n, h] ~ normal(mu_t_fall, sigma_t);
    }
  }
}

// Simulate a full observation from each value of the parameters
generated quantities {
  real obs_t_fall_post_pred[N_drops, N_heights];

  for (n in 1:N_drops) {
    for (h in 1:N_heights) {
      real mu_t_fall = t_fall(g, heights[h]);
      obs_t_fall_post_pred[n, h] = normal_rng(mu_t_fall, sigma_t);
    }
  }
}
