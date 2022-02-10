functions {
  vector tail_condition(vector y, vector theta, real[] x_r, int[] x_i) {
    vector[2] diff;

    if (!is_nan(y[1]) && !is_nan(y[2])) {
      diff[1] = inv_gamma_cdf(theta[1], exp(y[1]), exp(y[2])) - 0.01;
      diff[2] = inv_gamma_cdf(theta[2], exp(y[1]), exp(y[2])) - 0.99;
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
  real alpha_guess = square(delta * (theta_l + theta_u) / (theta_u - theta_l)) + 2;
  real beta_guess =  (alpha_guess - 1) * 0.5 * (theta_l + theta_u);
  vector[2] y_guess = [log(alpha_guess), log(beta_guess)]';
  vector[2] theta = [theta_l, theta_u]';
  real x_r[0];
  int x_i[0];
  
  print("Guesses:")
  print("  alpha = ", alpha_guess);
  print("  beta = ", beta_guess);
}

generated quantities {
  real alpha;
  real beta;
  {
    vector[2] y = algebra_solver(tail_condition, y_guess, theta, x_r, x_i);
    alpha = exp(y[1]);
    beta = exp(y[2]);
  }
}
