parameters {
  real x;
  real log_tau;
}

model {
  x ~ normal(0, exp(log_tau));
  log_tau ~ normal(0, 1);
}
