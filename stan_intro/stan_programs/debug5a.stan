parameters {
  real x[9];
}

model {
  for (n in 1:9)
    x ~ normal(0, 1);
}
