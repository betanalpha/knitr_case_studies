transformed data {
  row_vector a[3] = {1.0, 2.0, 3.0};
  vector b[3] = {1.0, 2.0, 3.0};
}

parameters {
  real x;
}

model {
  x ~ normal(a * b, 1);
}
