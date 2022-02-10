functions {
  vector tail_condition(vector y, vector theta, real[] x_r, int[] x_i) {
    real alpha = exp(y[1]);
    real delta = theta[2] - theta[1];
    vector[1] diff = [0.01 - 1 / (alpha * delta) * (log(2) - log1p(exp(-alpha * delta)))]';
    return diff;
  }
}

data {
  real theta_l;
  real theta_u;
}

transformed data {
  real alpha_guess = log(2) / (0.01 * (theta_u - theta_l));
  vector[1] y_guess = [log(alpha_guess)]';
  vector[2] theta = [theta_l, theta_u]';
  real x_r[0];
  int x_i[0];
  
  print("Guess:")
  print("  alpha = ", alpha_guess);
}

generated quantities {
  real alpha;
  {
    vector[1] y = algebra_solver(tail_condition, y_guess, theta, x_r, x_i);
    alpha = exp(y[1]);
  }
}
