parameters {
  vector<lower=0>[50] x;
}

model {
  x ~ cauchy(0, 1);
}
