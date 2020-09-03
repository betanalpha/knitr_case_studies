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

  y_obs ~ multi_normal_cholesky(rep_vector(0, N_obs), L_cov);
}
