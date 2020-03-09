############################################################
# Initial Setup
############################################################

library(rstan)
rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores()) # Parallelize chains

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

############################################################
# Bayesian Inference
############################################################

# Let's put all of this into the context of Bayesian inference!

# We first simulate an observation and save it to a file
N <- 1
simu_data <- list("N" = N)

simu <- stan(file='simulate_data.stan', data=simu_data,
             iter=1, chains=1, seed=4838282, algorithm="Fixed_param")

y <- array(extract(simu)$y[1,])
stan_rdump(c("N", "y"), file="simulation.data.R")

# Now we can read that data back in and use Hamiltonian
# Monte Carlo to estimate posterior expectation values
input_data <- read_rdump("simulation.data.R")
fit <- stan(file='fit_data.stan', data=input_data, seed=4938483)

# Check diagnostics
util$check_all_diagnostics(fit)

# That doesn't look good.  Let's investigate the samples
# to see what's going on.

partition <- util$partition_div(fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

plot(nondiv_params$mu, log(nondiv_params$sigma),
     col=c_dark_trans, pch=16, cex=0.8, 
     xlab="mu", ylab="log(sigma)")
points(div_params$mu, log(div_params$sigma),
       col=c_green_trans, pch=16, cex=0.8)

# What happens when we add more observations?
N <- 5
simu_data <- list("N" = N)

simu <- stan(file='simulate_data.stan', data=simu_data,
             iter=1, chains=1, seed=4838282, algorithm="Fixed_param")

y <- extract(simu)$y[1,]
stan_rdump(c("N", "y"), file="simulation.data.R")

# Now we can read that data back in and use Hamiltonian
# Monte Carlo to estimate posterior expectation values
input_data <- read_rdump("simulation.data.R")
fit <- stan(file='fit_data.stan', data=input_data, seed=4938483)

# Check diagnostics
util$check_all_diagnostics(fit)

# Much better
params = extract(fit)

plot(params$mu, log(params$sigma),
     col=c_dark_trans, pch=16, cex=0.8, 
     xlab="mu", ylab="log(sigma)")
