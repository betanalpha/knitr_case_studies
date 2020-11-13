functions {
  vector group_means(vector y, int N, int[] indiv_idxs, int N_indiv) {
    vector[N_indiv] counts = rep_vector(0, N_indiv);
    vector[N_indiv] means = rep_vector(0, N_indiv);
    for (n in 1:N) {
      int idx = indiv_idxs[n];
      real delta = y[n] - means[idx];
      counts[idx] += 1;
      means[idx] += delta / counts[idx];
    }
    return means;
  }
  
  real F_stat(vector y, int N, int[] indiv_idxs, int N_indiv) {
    real global_mean = 0;
    vector[N_indiv] counts = rep_vector(0, N_indiv);
    vector[N_indiv] means = rep_vector(0, N_indiv);
    vector[N_indiv] sum_squares = rep_vector(0, N_indiv);
    
    real between = 0;
    real within = 0;

    for (n in 1:N) {
      int idx = indiv_idxs[n];
      real delta;

      delta = y[n] - global_mean;
      global_mean += delta / n;

      counts[idx] += 1;
      delta = y[n] - means[idx];
      means[idx] += delta / counts[idx];
      sum_squares[idx] += (y[n] - means[idx]) * delta;
    }

    for (l in 1:N_indiv)
      between += counts[l] * square(means[l] - global_mean) / (N_indiv - 1);
    within = mean(sum_squares) / (N - N_indiv);

    return between / within;
  }
}

data {
  int<lower=0> N;      // Number of observations
  int y[N];            // Observed counts
  
  // Number of individual experiments
  int<lower=0> N_exp;
  // Experiment from which each observation is generated
  int<lower=1, upper=N_exp> exp_idx[N]; 
  
  int<lower=0, upper=N_exp> K_ncp;          // Number of noncentered individuals
  int<lower=1, upper=N_exp> ncp_idx[K_ncp]; // Index of noncentered individuals
  
  int<lower=0, upper=N_exp> K_cp;           // Number of centered individuals
  int<lower=1, upper=N_exp> cp_idx[K_cp];   // Index of noncentered individuals
  
  vector[N_exp] log_rho_I; // Nominal input for each experiment

  vector[N_exp] log_rho_T; // Temperature for each experiment
}

parameters {
  real mu_log_psi0;           // Analog-to-digital coefficients population location
  real<lower=0> tau_log_psi0; // Analog-to-digital coefficients population scale
  vector[K_ncp] log_psi0_ncp; // Non-centered individual analog-to-digital coefficients
  vector[K_cp]  log_psi0_cp;  // Ccentered individual analog-to-digital coefficients
  
  real gamma;                     // Exponential temperature dependence
}

transformed parameters {
  // Recentered individual parameters
  vector[N_exp] log_psi0;
  log_psi0[ncp_idx] = mu_log_psi0 + tau_log_psi0 *  log_psi0_ncp;
  log_psi0[cp_idx] =  log_psi0_cp;
}

model {
  vector[N_exp] log_psi = log_psi0 + gamma * log_rho_T;
  
  mu_log_psi0 ~ normal(0, log(2));
  tau_log_psi0 ~ normal(0, 0.5 * log(2));
  
  log_psi0_ncp ~ normal(0, 1);                     // Non-centered hierarchical model
  log_psi0_cp ~ normal(mu_log_psi0, tau_log_psi0); // Centered hierarchical model
  
  gamma ~ normal(1.25, 0.25);
  
  y ~ poisson_log(log_psi[exp_idx] + log_rho_I[exp_idx]); // Observational model
}

generated quantities {
  int y_post_pred[N] = poisson_log_rng(  log_psi0[exp_idx] 
                                       + gamma * log_rho_T[exp_idx] 
                                       + log_rho_I[exp_idx]);
  vector[N_exp] norm_ave = group_means(to_vector(y_post_pred) ./ exp(log_rho_I[exp_idx]), 
                                       N, exp_idx, N_exp);
  real F_exp = F_stat(to_vector(y_post_pred) ./ exp(log_rho_I[exp_idx]), 
                      N, exp_idx, N_exp);
}
