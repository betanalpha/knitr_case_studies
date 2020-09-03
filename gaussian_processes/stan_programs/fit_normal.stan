data {
  int<lower=1> N_predict;
  real x_predict[N_predict];

  int<lower=1> N_obs;
  real y_obs[N_obs];
  int<lower=1, upper=N_predict> observed_idx[N_obs];

  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

transformed data {
  matrix[N_predict, N_predict] cov =   cov_exp_quad(x_predict, alpha, rho)
                                     + diag_matrix(rep_vector(1e-10, N_predict));
  matrix[N_predict, N_predict] L_cov = cholesky_decompose(cov);
}

parameters {
  vector[N_predict] f_predict;
}

model {
  f_predict ~ multi_normal_cholesky(rep_vector(0, N_predict), L_cov);
  y_obs ~ normal(f_predict[observed_idx], sigma);
}

generated quantities {
  real y_predict[N_predict] = normal_rng(f_predict, sigma);
}
