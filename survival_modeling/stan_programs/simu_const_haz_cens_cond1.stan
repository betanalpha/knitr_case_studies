functions {
  // Inverse survival function
  real inv_survival(real u, real gamma) {
    return - log(u) / gamma;
  }
}

data {
  int<lower=1> N2;
  real<lower=0> gamma;
  
  real t1; // Left-censoring threshold
  real t2; // Interval-censoring left threshold
  real t3; // Interval-censoring right threshold/Right-censoring threshold
}

generated quantities {
  int<lower=0> N1 = 0;
  real<lower=0> obs_times[N2] = rep_array(0.0, N2);
  int<lower=0> N3 = 0;
  int<lower=0> N4 = 0;
  {
    int n = 1;
    while (1) {
      real u = uniform_rng(0, 1);
      real time = inv_survival(u, gamma);
      
      if (time < t1) {
        N1 += 1;
      } else if (time < t2) {
        obs_times[n] = time;
        n += 1;
      } else if (time < t3) {
        N3 += 1;
      } else {
        N4 += 1;
      }
     
      if (n > N2) break;  
    } 
  }
}
