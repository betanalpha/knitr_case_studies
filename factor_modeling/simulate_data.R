library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

util <- new.env()
source('stan_utility.R', local=util)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_trans <- c("#DCBCBC80")
c_dark_trans <- c("#8F272780")
c_green_trans <- c("#00FF0080")

par(family="CMU Serif", las=1, bty="l", cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 5))

################################################################################
# Simulate Single Factor Data
################################################################################

# Lots of data where monolithic centered parameterization is optimal
fit <- stan(file='stan_programs/simu_data.stan',
            iter=1, warmup=0, chains=1, seed=2838281, algorithm="Fixed_param")

simu_samples <- extract(fit)
sigma <- simu_samples$sigma[1]
N <- simu_samples$N[1]
K <- simu_samples$K[1]
obs_to_context <- simu_samples$obs_to_context[1,]
N_factor_levels <- simu_samples$N_factor_levels[1]
context_to_level <- simu_samples$context_to_level[1,]
y <- simu_samples$y[1,]

data <- list("sigma" = sigma,
             "N" = N, "K" = K, "obs_to_context" = obs_to_context,
             "N_factor_levels" = N_factor_levels, 
             "context_to_level" = context_to_level,
             "y" = y)

stan_rdump(names(data), file="data/single_factor1.data.R")

# More realistic data set where mixed parameterization is needed
fit <- stan(file='stan_programs/simu_data2.stan',
            iter=1, warmup=0, chains=1, seed=2838281, algorithm="Fixed_param")

simu_samples <- extract(fit)
sigma <- simu_samples$sigma[1]
N <- simu_samples$N[1]
K <- simu_samples$K[1]
obs_to_context <- simu_samples$obs_to_context[1,]
N_factor_levels <- simu_samples$N_factor_levels[1]
context_to_level <- simu_samples$context_to_level[1,]
y <- simu_samples$y[1,]

data <- list("sigma" = sigma,
             "N" = N, "K" = K, "obs_to_context" = obs_to_context,
             "N_factor_levels" = N_factor_levels, 
             "context_to_level" = context_to_level,
             "y" = y)

stan_rdump(names(data), file="data/single_factor2.data.R")

################################################################################
# Simulate Nested Factor Data
################################################################################

# Nested factors -- same data and model but different names
fit <- stan(file='stan_programs/simu_nested_data.stan',
            iter=1, warmup=0, chains=1, seed=8593938, algorithm="Fixed_param")

simu_samples <- extract(fit)

N <- simu_samples$N[1]
N_factor1_levels <- simu_samples$N_factor1_levels[1]
N_factor2_levels <- simu_samples$N_factor2_levels[1]
N_factor3_levels <- simu_samples$N_factor3_levels[1]

obs_to_factor3 <- simu_samples$obs_to_factor3[1,]
factor3_to_factor2 <- simu_samples$factor3_to_factor2[1,]
factor2_to_factor1 <- simu_samples$factor2_to_factor1[1,]

y <- simu_samples$y[1,]
sigma <- simu_samples$sigma[1]

names <- c("N", "N_factor1_levels", "N_factor2_levels", "N_factor3_levels",
           "obs_to_factor3", "factor3_to_factor2", "factor2_to_factor1",
           "y", "sigma")

stan_rdump(names, file="data/nested_multifactor.data.R")

################################################################################
# Simulate Overlapping Factor Data
################################################################################

simu <- stan(file='stan_programs/simu_overlapping_data.stan', seed=93585393,
             iter=1, warmup=0, chains=1, algorithm="Fixed_param")

simu_output <- extract(simu)

N <- simu_output$N[1]
N_main_factors <- simu_output$N_main_factors[1]
N_main_factor_levels <- simu_output$N_main_factor_levels[1,]
main_factor_level_idx <- simu_output$main_factor_level_idx[1,,]

y <- simu_output$y[1,]
sigma <- simu_output$sigma[1]

names <- c("N", "N_main_factors", "N_main_factor_levels", "main_factor_level_idx",
           "y", "sigma")

stan_rdump(names, file="data/overlapping_multifactor.data.R")

N <-348

main_factor_level_idx <- matrix(0, nrow=N_main_factors, N)
main_factor_level_idx[1,] <- sample(1:N_main_factor_levels[1], N, replace=TRUE)
main_factor_level_idx[2,] <- sample(1:N_main_factor_levels[2], N, replace=TRUE)
main_factor_level_idx[3,] <- rbinom(N, N_main_factor_levels[3] - 1, 0.8) + 1

names <- c("N", "N_main_factors", "N_main_factor_levels", "main_factor_level_idx")

stan_rdump(names, file="data/overlapping_multifactor_pred.data.R")

################################################################################
# Simulate Batch Variation Data
################################################################################

# Nested factors -- same data and model but different names
fit <- stan(file='stan_programs/simu_batch_data.stan',
            iter=1, warmup=0, chains=1, seed=8593938, algorithm="Fixed_param")

simu_samples <- extract(fit)

N <- simu_samples$N[1]
N_factor1_levels <- simu_samples$N_factor1_levels[1]
N_factor2_levels <- simu_samples$N_factor2_levels[1]
N_factor3_levels <- simu_samples$N_factor3_levels[1]

obs_to_factor3 <- simu_samples$obs_to_factor3[1,]
factor3_to_factor2 <- simu_samples$factor3_to_factor2[1,]
factor2_to_factor1 <- simu_samples$factor2_to_factor1[1,]

log_rho_I <- simu_samples$log_rho_I[1,]
y <- simu_samples$y[1,]

names <- c("N", "N_factor1_levels", "N_factor2_levels", "N_factor3_levels",
           "obs_to_factor3", "factor3_to_factor2", "factor2_to_factor1",
           "y", "log_rho_I")

stan_rdump(names, file="data/batch.data.R")

################################################################################
# Simulate Quality Assurance Data
################################################################################

simu <- stan(file='stan_programs/simu_qa_data.stan', seed=93585392,
             iter=1, warmup=0, chains=1, algorithm="Fixed_param")

simu_output <- extract(simu)

N <- simu_output$N[1]
N_factor_configs <- simu_output$N_factor_configs[1]
N_samples <- simu_output$N_samples[1]
factor_config_idx <- simu_output$factor_config_idx[1,]
y <- simu_output$y[1,]
t <- simu_output$t[1,]

N_location <- simu_output$N_main_factor_levels[1,1]
location_idx <- simu_output$main_factor_level_idx[1,1,]

N_method <- simu_output$N_main_factor_levels[1,2]
method_idx <- simu_output$main_factor_level_idx[1,2,]

N_source <- simu_output$N_main_factor_levels[1,3]
source_idx <- simu_output$main_factor_level_idx[1,3,]

names <- c("N", "N_samples", "N_factor_configs",
           "factor_config_idx", "y", "t",
           "N_location", "location_idx",
           "N_method", "method_idx",
           "N_source", "source_idx")

stan_rdump(names, file="data/qa.data.R")
