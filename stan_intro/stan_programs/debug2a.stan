transformed data {
  real a = 5;
}

parameters {
  real x;
}

model {
  real b;
  x ~ normal(a * b, 1);
}
