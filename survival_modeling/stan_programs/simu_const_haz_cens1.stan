functions {
  // Inverse survival function
  real inv_survival(real u, real gamma) {
    return - log(u) / gamma;
  }
}

data {
  int<lower=1> N;
  real<lower=0> gamma;
  
  real t1; // Left-censoring threshold
  real t2; // Interval-censoring left threshold
  real t3; // Interval-censoring right threshold/Right-censoring threshold
}

transformed data {
  real<lower=0> all_obs_times[N] = rep_array(0.0, N);
  int<lower=0> M1 = 0;
  int<lower=0> M2 = 0;
  int<lower=0> M3 = 0;
  int<lower=0> M4 = 0;
  
  for (n in 1:N) {
    real u = uniform_rng(0, 1);
    real time = inv_survival(u, gamma);
      
    if (time < t1) {
      M1 += 1;
    } else if (time < t2) {
      M2 += 1;
      all_obs_times[M2] = time;
    } else if (time < t3) {
      M3 += 1;
    } else {
      M4 += 1;
    }
  }  
}

generated quantities {
  int<lower=0> N1 = M1;
    
  int<lower=0> N2 = M2;
  real<lower=0> obs_times[M2] = all_obs_times[1:M2];
  
  int<lower=0> N3 = M3;
  int<lower=0> N4 = M4;
}
