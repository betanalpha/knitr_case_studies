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

  // Nominal input for each experiment
  vector[N_exp] log_rho_I;
}

generated quantities {
  vector[N_exp] norm_ave = group_means(to_vector(y) ./ exp(log_rho_I[exp_idx]), N, exp_idx, N_exp);
  real F_exp = F_stat(to_vector(y) ./ exp(log_rho_I[exp_idx]), N, exp_idx, N_exp);
}
