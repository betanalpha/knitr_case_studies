parameters {
  vector[50] x;
}

model {
  x ~ cauchy(0, 1);
}
