transformed data {
  real a = 5;
}

parameters {
  real x;
}

model {
  real b = 4;
  x ~ normal(a * b, 1);
}
