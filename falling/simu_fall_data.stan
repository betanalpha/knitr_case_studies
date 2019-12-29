functions {
  // Fall time in vacuum
  real t_fall(real g, real h, real v0) {
    return (v0 + sqrt(square(v0) + 2 * g * h)) / g;
  }
  
  // Fall time in fluid
  real t_fall_air(real g, real h, real v0, real k, real m) {
    real vT = sqrt(m * g / k);
    
    if (v0 > 0) {
      // Positive initial velocity
      real a1 = acosh(  exp(  g * h / square(vT)
                            + 0.5 * log(1 + square(v0 / vT)) ));
      real a2 = atan(v0 / vT);
      return (vT / g) * (a1 + a2);
    } else if (v0 > -vT) {
      // Negative but subterminal initial velocity
      real a1 = acosh(  exp( g * h / square(vT)
                            - 0.5 * log(1 - square(v0 / vT)) ));
      real a2 = atanh(v0 / vT);
      return (vT / g) * (a1 + a2);
    } else {
      // Negative and superterminal initial velocity
      real a1 = asinh(  exp(  g * h / square(vT)
                             -0.5 * log(square(v0 / vT) - 1) )); 
      real a2 = atanh(vT / v0);
      return (vT / g) * (a1 + a2);
    }
  }
}

data {
  int N_heights;
  real<lower=0> heights[N_heights]; // meters
  int N_drops;
  real m; // kilograms
}

transformed data {
  real sigma_t = 0.02; // seconds
  real sigma_v = 0.5;  // meters per second
  real g = 9.806;      // meters per second squared
  real k = 0.0012;     // kilograms per meter
}

generated quantities {
  real obs_t_fall[N_drops, N_heights];
  real v0[N_drops];

  for (n in 1:N_drops) {
    v0[n] = normal_rng(0, sigma_v);
    for (h in 1:N_heights) {
      real mu_t_fall_air = t_fall_air(g, heights[h], v0[n], k, m);
      obs_t_fall[n, h] = normal_rng(mu_t_fall_air, sigma_t);
    }
  }
}
