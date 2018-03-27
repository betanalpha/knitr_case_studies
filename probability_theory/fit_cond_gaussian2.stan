transformed data {
  real x1 = 1;

  vector[2] mu = [2, -1]';

  real sigma1 = 0.75;
  real sigma2 = 1.2;
  real rho = 0.7;

  matrix[2,2] Sigma;
  Sigma[1][1] = sigma1 * sigma1;
  Sigma[1][2] = rho * sigma1 * sigma2;
  Sigma[2][1] = rho * sigma1 * sigma2;
  Sigma[2][2] = sigma2 * sigma2;
}

parameters {
  real x2;
}

model {
  vector[2] x = [x1, x2]';
  x ~ multi_normal(mu, Sigma);
}
