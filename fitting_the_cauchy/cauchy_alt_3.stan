parameters {
  vector<lower=0, upper=1>[50] x_tilde;
}

transformed parameters {
vector[50] x = tan(pi() * (x_tilde - 0.5));
}

model {
  // Implicit uniform prior on x_tilde
}

generated quantities {
  real I = fabs(x[1]) < 1 ? 1 : 0;
}
