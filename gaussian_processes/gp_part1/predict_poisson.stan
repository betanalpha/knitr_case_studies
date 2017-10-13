data {
  int<lower=1> N_predict;
  real x_predict[N_predict];

  int<lower=1> N_observed;
  int<lower=1, upper=N_predict> observed_idx[N_observed];
  int y_observed[N_observed];

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
  vector[N_predict] log_f_predict = L_cov * f_tilde;
}

model {
  f_tilde ~ normal(0, 1);
  y_observed ~ poisson_log(log_f_predict[observed_idx]);
}

generated quantities {
  vector[N_predict] f_predict = exp(log_f_predict);
  vector[N_predict] y_predict;
  for (n in 1:N_predict)
    y_predict[n] = poisson_log_rng(log_f_predict[n]);
}
