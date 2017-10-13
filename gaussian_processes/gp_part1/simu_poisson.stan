data {
  int<lower=1> N;
  real x[N];

  real<lower=0> rho;
  real<lower=0> alpha;
}

transformed data {
  matrix[N, N] cov =   cov_exp_quad(x, alpha, rho)
                     + diag_matrix(rep_vector(1e-10, N));
  matrix[N, N] L_cov = cholesky_decompose(cov);
}

parameters {}
model {}

generated quantities {
  vector[N] f = multi_normal_cholesky_rng(rep_vector(0, N), L_cov);
  int y[N];
  for (n in 1:N)
    y[n] = poisson_log_rng(f[n]);
}
