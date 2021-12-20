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
  // Baseline population exposure rate
  real<lower=0> gamma;
  
  // Susceptibilty of individuals employeed at first company
  // Units of inverse exposure
  real<lower=0> alpha1;
  
  // Susceptibilty of individuals employeed at second company
  // Units of inverse exposure
  real<lower=0> alpha2;
  
  // Treatment effect, units of inverse exposure
  real<lower=-alpha2> zeta;
}

model {
  gamma ~ normal(0, 1);
  
  alpha1 ~ normal(0, 1);
  alpha2 ~ normal(0, 1);
  zeta ~ normal(0, 1);
  
  x1 ~ exponential(gamma * alpha1);
  y1 ~ bernoulli(inf_probs(alpha1, x1));
  
  x2 ~ exponential(gamma * alpha2);
  y2 ~ bernoulli(inf_probs(alpha2 + zeta, x2));
}

generated quantities {
  real x1_pred[N1]; // Covariate retrodiction for first company
  real y1_pred[N1]; // Conditional variate retrodiction for first company
  real x2_pred[N2]; // Covariate retrodiction for second company
  real y2_pred[N2]; // Conditional variate retrodiction for second company
  
  real p_nt_grid[N_grid]; // Inferred probability verses exposure with no treatment
  real p_t_grid[N_grid];  // Inferred probability verses exposure with treatment

  for (n in 1:N1) {
    x1_pred[n] = exponential_rng(gamma * alpha1);
    y1_pred[n] = bernoulli_rng(inf_prob(alpha1, x1[n]));
  }
  
  for (n in 1:N2) {
    x2_pred[n] = exponential_rng(gamma * alpha2);
    y2_pred[n] = bernoulli_rng(inf_prob(alpha2 + zeta, x2[n]));
  }

  for (n in 1:N_grid) {
    p_nt_grid[n] = inf_prob(alpha1, x_grid[n]);
    p_t_grid[n] = inf_prob(alpha2 + zeta, x_grid[n]);
  }
}
