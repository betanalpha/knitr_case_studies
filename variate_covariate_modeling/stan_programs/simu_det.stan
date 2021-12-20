functions {
  // Location-dispersion paramaeterization of gamma family
  real gamma_ld_lpdf(vector y, vector mu, real phi) {
    real beta = inv(phi);
    vector[num_elements(mu)] alpha = beta * mu;
    return gamma_lpdf(y | alpha, beta);
  }
  
  // Location-dispersion paramaeterization of gamma family
  real gamma_ld_rng(real mu, real phi) {
    real beta = inv(phi);
    real alpha = beta * mu;
    return gamma_rng(alpha, beta);
  }
}

data {
  int<lower=1> N1;
  int<lower=1> N2;
}

transformed data {
  real lambda = 0.1;
  real psi = 0.25;
  
  real mu1 = 30;
  real phi1 = 15;
  
  real mu2 = 70;
  real phi2 = 1;
}

generated quantities {
  real x1[N1];
  real y1[N1];
  
  real x2[N2];
  real y2[N2];

  for (n in 1:N1) {
    x1[n] = gamma_ld_rng(mu1, phi1);
    y1[n] = gamma_ld_rng(lambda * x1[n], psi);
  }
  
  for (n in 1:N2) {
    x2[n] = gamma_ld_rng(mu2, phi2);
    y2[n] = gamma_ld_rng(lambda * x2[n], psi);
  }
}
