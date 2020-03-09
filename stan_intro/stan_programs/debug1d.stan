transformed data {
  row_vector[3] a = [1.0, 2.0, 3.0];
  vector[3] b = [1.0, 2.0, 3.0];
}

parameters {
  real x;
}

model {
  x ~ normal(a * b, 1);
}
