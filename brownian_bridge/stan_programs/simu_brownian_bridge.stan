functions {
  real bridge_rng(real t, real t_left, real x_left, real t_right, real x_right) {
    real mu = x_left + (t - t_left) / (t_right - t_left) * (x_right - x_left);
    real sigma2 = (t_right - t) * (t - t_left) / (t_right - t_left);
    return normal_rng(mu, sqrt(sigma2));
  }
}

data {
  int<lower=1> N;
  real ts[N];
  
  real x_left;
  real x_right;
  
  int left_idxs[N - 2];
  int center_idxs[N - 2];
  int right_idxs[N - 2];
}

generated quantities {
  real xs[N];
  
  xs[1] = x_left;
  xs[N] = x_right;
  
  for (n in 1:(N - 2)) {
    xs[center_idxs[n]] = bridge_rng(ts[center_idxs[n]],
                                    ts[left_idxs[n]],  xs[left_idxs[n]],
                                    ts[right_idxs[n]], xs[right_idxs[n]]);
  }
  
}
