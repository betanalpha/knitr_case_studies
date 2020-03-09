transformed data {
  real a[3] = {1.0, 2.0, 3.0};
  real b[3] = {1.0, 2.0, 3.0};
}

parameters {
  real x;
}

model {
  x ~ normal(a * b, 1);
}
