data {
  int<lower=1> N1;
  int<lower=1> N2;
}

transformed data {
  real<lower=0> A = 3;
  real phi = 4;
  
  real<lower=0> omega1 = 2;
  real<lower=0> alpha1 = 3;
  real<lower=0> beta1 = 0.25;
  
  real<lower=0> omega2 = 1;
  real<lower=0> alpha2 = 38;
  real<lower=0> beta2 = 4.5;
  
  real gamma = 0.5;
  real sigma = 0.5;
}

generated quantities {
  real x1[N1];
  real y1[N1];
  
  real x2[N2];
  real y2[N2];

  for (n in 1:N1) {
    x1[n] = gamma_rng(alpha1, beta1 + gamma * omega1);
    y1[n] = normal_rng(A * sin(omega1 * x1[n] + phi), sigma);
  }
  
  for (n in 1:N2) {
    x2[n] = gamma_rng(alpha2, beta2 + gamma * omega2);
    y2[n] = normal_rng(A * sin(omega2 * x2[n] + phi), sigma);
  }
}
