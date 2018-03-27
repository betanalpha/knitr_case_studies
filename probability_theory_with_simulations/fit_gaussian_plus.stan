parameters {
  real x;
}

model {
  x ~ normal(0, 1);
}

generated quantities {
  real sq_dev = square(x - 0); // True mean is zero
  real P = -2 < x && x < 1 ? 1 : 0;
}
