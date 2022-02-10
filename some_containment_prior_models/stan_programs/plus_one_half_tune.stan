functions {
  real cdf(real x, real mu, real sigma) {
    real denom = sqrt(sigma * x);
    return Phi((x - mu) / denom) - exp(2 * mu / sigma) * Phi(- (x + mu) / denom);
  }
    
  vector tail_condition(vector y, vector theta, real[] x_r, int[] x_i) {
    vector[2] diff;
    if (!is_nan(y[1]) && !is_nan(y[2])) {
      diff[1] = cdf(theta[1], y[1], exp(y[2])) - 0.01;
      diff[2] = cdf(theta[2], y[1], exp(y[2])) - 0.99;
    } else {
      diff[1] = 0;
      diff[2] = 0;
    }
    return diff;
  }
}

data {
  real theta_l;
  real theta_u;
  real delta;
}

transformed data {
  real sigma_guess = 
      0.25 * (theta_l + theta_u) 
    * (sqrt(1 + 4 * square( (theta_u - theta_l) / (delta * (theta_l + theta_u)) ) ) - 1);
  real mu_guess = 0.5 * (theta_l + theta_u) - sigma_guess;
  vector[2] y_guess = [mu_guess, log(sigma_guess)]';
  vector[2] theta = [theta_l, theta_u]';
  real x_r[0];
  int x_i[0];
  
  print("Guesses:")
  print("  mu = ", mu_guess);
  print("  sigma = ", sigma_guess);
}

generated quantities {
  real mu;
  real sigma;
  {
    vector[2] y = algebra_solver(tail_condition, y_guess, theta, x_r, x_i);
    mu = y[1];
    sigma = exp(y[2]);
  }
}
