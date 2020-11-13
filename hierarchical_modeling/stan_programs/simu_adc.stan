transformed data {
  real mu_log_psi0 = 0.5;
  real<lower=0> tau_log_psi0 = 0.5;
  real gamma = 1.25;
  int N_prime = 1100;
  int N_exp_prime = 21;
}

generated quantities {
  int<lower=0> N = N_prime;
  int<lower=0> N_exp = N_exp_prime;

  real log_rho_I[N_exp_prime] = {0.00, 0.1, 0.30, 0.45, 0.60, 0.75, 0.90, 
                                 1.05, 1.20, 1.35, 1.50, 1.65, 1.80, 1.95, 
                                 2.10, 2.25, 2.40, 2.55, 2.70, 2.85, 3.00};
  real log_rho_T[N_exp_prime];
  
  int<lower=1, upper=N_exp> exp_idx[N_prime]; 
  int y[N_prime];
  
  {
    real log_psi0[N_exp_prime];
    vector[N_exp_prime] exp_probs = [0.090090090, 0.009009009, 0.009009009, 0.009009009, 
                                     0.090090090, 0.090090090, 0.009009009, 0.009009009, 
                                     0.009009009, 0.009009009, 0.090090090, 0.090090090, 
                                     0.090090090, 0.009009009, 0.090090090, 0.009009009, 
                                     0.090090090, 0.090090090, 0.009009009, 0.090090090, 
                                     0.009009009]';

    for (k in 1:N_exp) {
      log_psi0[k] = normal_rng(mu_log_psi0, tau_log_psi0);
      log_rho_T[k] = normal_rng(0, log(2));
    }

    for (n in 1:N) {
      exp_idx[n] = categorical_rng(exp_probs);
      y[n] = poisson_log_rng(  log_psi0[exp_idx[n]] 
                             + gamma * log_rho_T[exp_idx[n]]
                             + log_rho_I[exp_idx[n]]);
    }

  }
}
