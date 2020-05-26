transformed data {
 real theta = pi() / 4;
 real phi = 0;
}

parameters {
  real x;
  real y;
  real z;
}

model {
  target += - 0.5 * square( (y - x * x - x + 3) / 0.25)
            - 0.5 * square(x / 2) - 0.5 * square(y / 2) - 0.5 * square(z / 2);
}

generated quantities {
  real x_rot = cos(phi) * (cos(theta) * x + sin(theta) * z) + sin(phi) * y;
  real y_rot = -sin(phi) * (cos(theta) * x + sin(theta) * z) + cos(phi) * y;
  real z_rot = - sin(theta) * x + cos(theta) * z;
}
