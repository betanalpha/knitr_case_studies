functions {
  // Infection probability as a function of susceptibilty and exposure
  real inf_prob(real alpha, real x) {
    return -expm1(-alpha * x); // 1 - exp(-alpha * x)
  }
  
  vector inf_probs(real alpha, vector x) {
    int N = num_elements(x);
    vector[N] p;
    for (n in 1:N) p[n] = inf_prob(alpha, x[n]);
    return p;
  }
}

data {
  int<lower=1> N1; // Number of observations at first comapny
  vector[N1] x1;   // Covariate observations at first company
  int y1[N1];      // Variate observations at first company
  
  int<lower=1> N2; // Number of observations at second comapny
  vector[N2] x2;   // Covariate observations at second company
  int y2[N2];      // Variate observations at second company
  
  int<lower=1> N_grid; // Number of grid points for quantifying functional behavior
  vector[N_grid] x_grid; // Grid points for quantifying functional behavior
}

parameters { 
  real<lower=0> alpha;     // Susceptibilty, units of inverse exposure
  real<lower=-alpha> zeta; // Treatment effect, units of inverse exposure
}

model {
  alpha ~ normal(0, 1);
  zeta ~ normal(0, 1);
  
  y1 ~ bernoulli(inf_probs(alpha, x1));
  y2 ~ bernoulli(inf_probs(alpha + zeta, x2));
}

generated quantities {
  real y1_pred[N1]; // Conditional variate retrodiction for first company
  real y2_pred[N2]; // Conditional variate retrodiction for second company
  
  real p_nt_grid[N_grid]; // Inferred probability verses exposure with no treatment
  real p_t_grid[N_grid];  // Inferred probability verses exposure with treatment

  for (n in 1:N1) {
    y1_pred[n] = bernoulli_rng(inf_prob(alpha, x1[n]));
  }
  
  for (n in 1:N2) {
    y2_pred[n] = bernoulli_rng(inf_prob(alpha + zeta, x2[n]));
  }

  for (n in 1:N_grid) {
    p_nt_grid[n] = inf_prob(alpha, x_grid[n]);
    p_t_grid[n] = inf_prob(alpha + zeta, x_grid[n]);
  }
}
