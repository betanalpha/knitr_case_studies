transformed data {
  real a = 5;
}

parameters {
  real x;
}

model {
  real b;
  print("a * b = ", a * b, ", a = ", a, ", b = ", b);
  x ~ normal(a * b, 1);
}
