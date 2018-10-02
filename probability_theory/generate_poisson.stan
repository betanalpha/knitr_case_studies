// Values provided externally
data {
  real l;
}

generated quantities {
  // Generate an exact sample from a distribution
  // specified by a Poisson probabilty mass function
  int x = poisson_rng(l);
}
