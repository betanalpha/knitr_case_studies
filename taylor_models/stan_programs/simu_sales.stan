transformed data {
  int<lower=0> N = 230; // Number of observations
}

generated quantities {
  real t[N]; // Temperature in F
  real p[N]; // Price in cents
  real b[N]; // Hotel bookings in counts
  
  real s[N]; // Sales

  for (n in 1:N) {
    t[n] = normal_rng(78, 4);
    p[n] = normal_rng(300, 15);
    b[n] = normal_rng(1000 + 50 * (t[n] - 75), 200);
    
    s[n] = normal_rng(  2910
                      + 0.75 * t[n] - 1.5 * square(t[n] - 82) - 35
                      - 5 * (p[n] - 300)
                      + 0.2 * (p[n] - 300) * (t[n] - 75)
                      + 0.33 * (b[n] - 1000), 100);
  }
}
