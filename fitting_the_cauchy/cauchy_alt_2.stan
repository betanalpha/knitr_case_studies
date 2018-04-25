parameters {
  vector[50] x_a;
  vector<lower=0>[50] x_b;
}

transformed parameters {
  vector[50] x = x_a .* sqrt(x_b);
}

model {
  x_a ~ normal(0, 1);
  x_b ~ inv_gamma(0.5, 0.5);
}

generated quantities {
  real I = fabs(x[1]) < 1 ? 1 : 0;
}
