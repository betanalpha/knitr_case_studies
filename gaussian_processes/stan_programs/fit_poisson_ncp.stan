data {
  int<lower=1> N_predict;
  real x_predict[N_predict];

  int<lower=1> N_obs;
  int y_obs[N_obs];
  int<lower=1, upper=N_predict> observed_idx[N_obs];

  real<lower=0> rho;
  real<lower=0> alpha;
}

transformed data {
  matrix[N_predict, N_predict] cov =   cov_exp_quad(x_predict, alpha, rho)
                                     + diag_matrix(rep_vector(1e-10, N_predict));
  matrix[N_predict, N_predict] L_cov = cholesky_decompose(cov);
}

parameters {
  vector[N_predict] f_tilde;
}

transformed parameters {
  vector[N_predict] f_predict = L_cov * f_tilde;
}

model {
  f_tilde ~ normal(0, 1);
  y_obs ~ poisson_log(f_predict[observed_idx]);
}

generated quantities {
  int y_predict[N_predict] = poisson_log_rng(f_predict);
}
