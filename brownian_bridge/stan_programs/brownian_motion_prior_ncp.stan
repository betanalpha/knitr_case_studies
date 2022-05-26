data {
  int<lower=1> N;
  real ts[N];
  
  real x_left;
}

parameters {
  real xs_forward_ncp[N - 1];
}

transformed parameters {
  real xs[N]; 

  xs[1] = x_left;

  for (n in 2:N) {  
    xs[n] 
      = xs[n - 1] + sqrt(ts[n] - ts[n - 1]) * xs_forward_ncp[n - 1];
  }
}

model {
  xs_forward_ncp ~ normal(0, 1);
}
