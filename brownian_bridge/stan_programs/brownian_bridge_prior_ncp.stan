functions {
  real bridge_lpdf(real x, real t, real t_left, real x_left, real t_right, real x_right) {
    real mu = x_left + (t - t_left) / (t_right - t_left) * (x_right - x_left);
    real sigma2 = (t_right - t) * (t - t_left) / (t_right - t_left);
    return normal_lpdf(x | mu, sqrt(sigma2));
  }
  
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

parameters {
  real interior_x_tildes[N - 2];
}

transformed parameters {
  real xs[N];
  
  xs[1] = x_left;
  xs[N] = x_right;
    
  for (n in 1:(N - 2)) {
    real t = ts[center_idxs[n]];
    real t_l = ts[left_idxs[n]];
    real x_l = xs[left_idxs[n]];
    real t_r = ts[right_idxs[n]];
    real x_r =  xs[right_idxs[n]];
    
    real mu = x_l + (t - t_l) / (t_r - t_l) * (x_r - x_l);
    real sigma2 = (t_r - t) * (t - t_l) / (t_r - t_l);
    
    xs[center_idxs[n]] = mu + sqrt(sigma2) * interior_x_tildes[n];
  }
}

model {
  interior_x_tildes ~ normal(0, 1);
}
