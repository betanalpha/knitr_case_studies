functions {
  // Differences between Beta tail probabilities and target probabilities
  vector tail_delta(vector y, vector theta, real[] x_r, int[] x_i) {
    vector[2] deltas;
    deltas[1] = beta_cdf(theta[1], exp(y[1]), exp(y[2])) - 0.01;
    deltas[2] = 1 - beta_cdf(theta[2], exp(y[1]), exp(y[2])) - 0.01;
    return deltas;
  }
}

transformed data {
  vector[2] y_guess = [log(5), log(5)]'; // Initial guess at Beta parameters
  vector[2] theta = [0.1, 0.9]';  // Target quantile
  vector[2] y;
  real x_r[0];
  int x_i[0];

  // Find Beta density parameters that ensure
  // 1% probabilty below 0.1 and 1% probabilty above 0.99
  y = algebra_solver(tail_delta, y_guess, theta, x_r, x_i);

  print("alpha = ", exp(y[1]));
  print("beta = ", exp(y[2]));
}

generated quantities {
  real alpha = exp(y[1]);
  real beta = exp(y[2]);
}
