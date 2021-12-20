functions {
  real inf_prob(real alpha, real x) {
    return -expm1(-alpha * x); // 1 - exp(-alpha * x)
  }
  
  vector inf_probs(real alpha, vector x) {
    int N = num_elements(x);
    vector[N] p;
    for (n in 1:N) p[n] = inf_prob(alpha, x[n]);
    return p;
  }
}

data {
  int<lower=1> N1;
  int<lower=1> N2;
}

transformed data {
  real alpha = 0.5;
  real zeta = -0.35;
  real phi1 = 0.25;
  real phi2 = 0.75;
}

generated quantities {
  real x1[N1];
  int y1[N1];
  
  real x2[N2];
  int y2[N2];

  for (n in 1:N1) {
    x1[n] = exponential_rng(phi1);
    y1[n] = bernoulli_rng(inf_prob(alpha, x1[n]));
  }
  
  for (n in 1:N2) {
    x2[n] = exponential_rng(phi2);
    y2[n] = bernoulli_rng(inf_prob(alpha, x2[n]));
  }
}
