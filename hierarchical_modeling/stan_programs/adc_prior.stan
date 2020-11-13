generated quantities {
  real<lower=0> psi;
  {
    real mu_log_psi0 = normal_rng(0, 0.5 * log(1.25));
    real tau_log_psi0 = fabs(normal_rng(0, log(1.25)));
    real log_psi0 = normal_rng(mu_log_psi0, tau_log_psi0);
    real log_rho_T = normal_rng(0, log(2));
    real gamma = normal_rng(1.25, 0.25);
    psi = exp(log_psi0 + gamma * log_rho_T);
  }
}
