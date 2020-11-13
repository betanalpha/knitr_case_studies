data {
  int<lower=0> N;      // Number of observations
  vector[N] y;         // Observations
  real<lower=0> sigma; // Measurement variability
  
  // Number of individual contexts in hierarchy
  int<lower=0> K;
  // Individual context from which each observation is generated
  int<lower=1, upper=K> indiv_idx[N]; 
  
  int<lower=0, upper=K> K_ncp;          // Number of noncentered individuals
  int<lower=1, upper=K> ncp_idx[K_ncp]; // Index of noncentered individuals
  
  int<lower=0, upper=K> K_cp;           // Number of centered individuals
  int<lower=1, upper=K> cp_idx[K_cp];   // Index of noncentered individuals
}

parameters {
  real mu;                  // Population location
  real<lower=0> tau;        // Population scale
  vector[K_ncp] theta_ncp;  // Non-centered individual parameters
  vector[K_cp]  theta_cp;   // Ccentered individual parameters
}

transformed parameters {
  // Recentered individual parameters
  vector[K] theta;
  theta[ncp_idx] = mu + tau * theta_ncp;
  theta[cp_idx] = theta_cp;
}

model {
  mu ~ normal(0, 5);          // Prior model
  tau ~ normal(0, 5);         // Prior model
  
  theta_ncp ~ normal(0, 1);   // Non-centered hierarchical model
  theta_cp ~ normal(mu, tau); // Centered hierarchical model
  
  y ~ normal(theta[indiv_idx], sigma); // Observational model
}
