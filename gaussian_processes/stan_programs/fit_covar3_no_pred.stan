data {
  int<lower=1> N_obs;
  real x_obs[N_obs];
  vector[N_obs] y_obs;
}

parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

model {
  matrix[N_obs, N_obs] cov =   cov_exp_quad(x_obs, alpha, rho)
                             + diag_matrix(rep_vector(square(sigma), N_obs));
  matrix[N_obs, N_obs] L_cov = cholesky_decompose(cov);

  rho ~ inv_gamma(4.6, 22.1); // P[rho < 2.0] \approx 0.01, P[rho > 20] \approx 0.01
  alpha ~ normal(0, 2);
  sigma ~ normal(0, 1);

  y_obs ~ multi_normal_cholesky(rep_vector(0, N_obs), L_cov);
}
