parameters {
  vector[50] x;
}

model {
  x ~ cauchy(0, 1);
}

generated quantities {
  real I = fabs(x[1]) < 1 ? 1 : 0;
}
