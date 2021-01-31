data {
  int<lower=0> N_fine_t;
  vector[N_fine_t] fine_t;
}

parameters {    
  // Zeroth-order contribution (Baseline)
  real baseline;
  
  // Logit probability at time zero
  real alpha;
}

model {
  baseline ~ normal(-3, 2);
  alpha ~ normal(4, 2);
}

generated quantities {
  vector[N_fine_t] p = inv_logit(alpha + baseline * fine_t / 365.0);
}
