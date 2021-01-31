transformed data {
  real factor_mu = 5;
  real factor_tau = 2;
  
  //int N_latent = 2025;
  int N_latent = 810;
  int K_latent = 81;
  int L_latent = 9;
  
  vector[K_latent] context_probs = rep_vector(1.0 / K_latent, K_latent);
  vector[L_latent] level_probs = rep_vector(1.0 / L_latent, L_latent);
  
  real factor_theta[L_latent] 
    = normal_rng(rep_vector(factor_mu, L_latent), factor_tau);

  real<lower=0> residual_tau[L_latent]
    = fabs(normal_rng(rep_vector(0, L_latent), 3));
}

generated quantities {
  real sigma = 5;
    
  int<lower=1> N = N_latent;

  // Number of individual contexts
  int<lower=0> K = K_latent;
  int<lower=1, upper=K> obs_to_context[N_latent];
  
  // Factor information
  int<lower=0> N_factor_levels = L_latent;
  int<lower=1, upper=N_factor_levels> context_to_level[K_latent];
  
  vector[N_latent] y;  
  
  { 
    vector[K_latent] context_theta;
    
    for (k in 1:K) {
      int l = categorical_rng(level_probs);
      context_to_level[k] = l;
      context_theta[k] = normal_rng(factor_theta[l], residual_tau[l]);
    }

    for (n in 1:N) {
      int k = categorical_rng(context_probs);
      obs_to_context[n] = k;
      y[n] = normal_rng(context_theta[k], sigma);
    }  
  }
}
