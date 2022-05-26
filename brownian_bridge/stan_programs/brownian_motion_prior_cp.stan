data {
  int<lower=1> N;
  real ts[N];
  
  real x_left;
}

parameters {
  real xs_forward[N - 1];
}

transformed parameters {
  real xs[N] = append_array(rep_array(x_left, 1), xs_forward);
}

model {
  for (n in 2:N) {
    xs[n] ~ normal(xs[n - 1], sqrt(ts[n] - ts[n - 1]) );
  }
}
