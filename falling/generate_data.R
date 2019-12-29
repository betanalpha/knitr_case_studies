############################################################
# Initial setup
############################################################

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

util <- new.env()
source('stan_utility.R', local=util)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_trans <- c("#DCBCBC20")
c_dark_trans <- c("#8F272720")

############################################################
# Simulate detector calibration data
############################################################

simu_data <- list("N_calibrations" = 500)

simu_fit <- stan(file="simu_calib_data.stan", data=simu_data, iter=1,
                 chains=1, seed=194838, algorithm="Fixed_param")

calib_data <- list("N_calibrations" = simu_data$N_calibrations,
                   "obs_t_calib" = extract(simu_fit)$obs_t[1,])

N_calibrations <- calib_data$N_calibrations
obs_t_calib <- calib_data$obs_t_calib

stan_rdump(c("N_calibrations", "obs_t_calib"), file="det_calib.data.R")

############################################################
# Simulate falling data with heavy ball
############################################################

simu_data <- list("N_heights" = 6, "heights" = c(1, 2, 3, 4, 5, 6), "N_drops" = 200)

simu_data[["m"]] <- 0.5

simu_fit <- stan(file="simu_fall_data.stan", data=simu_data, iter=1,
                 chains=1, seed=194838, algorithm="Fixed_param")

obs_t_fall <- extract(simu_fit)$obs_t_fall[1,,]

data <- simu_data
data[["obs_t_fall"]] <- obs_t_fall
data[["m"]] <- rep(simu_data$m, simu_data$N_drops)

N_heights <- data$N_heights
heights <- data$heights
N_drops <- data$N_drops
obs_t_fall <- data$obs_t_fall
m <- data$m

stan_rdump(c("N_heights", "heights", "N_drops", "obs_t_fall", "m"), file="falling.data.R")

# Find extreme initial velocities
v0s <- extract(simu_fit)$v0[1,]
which(abs(v0s) < 0.01)
which(v0s > 1.1)
which(v0s < -1)
# Neutral 191, Positive 132, Negative 83

############################################################
# Simulate falling data with light ball
############################################################

simu_data[["m"]] <- 0.002

simu_fit <- stan(file="simu_fall_data.stan", data=simu_data, iter=1,
                 chains=1, seed=588337, algorithm="Fixed_param")

obs_t_fall <- extract(simu_fit)$obs_t_fall[1,,]

light_data <- simu_data
light_data[["obs_t_fall"]] <- obs_t_fall
light_data[["m"]] <- rep(simu_data$m, simu_data$N_drops)

obs_t_fall <- light_data$obs_t_fall
m <- light_data$m

stan_rdump(c("N_heights", "heights", "N_drops", "obs_t_fall", "m"), file="falling_light.data.R")

# Find extreme initial velocities
v0s <- extract(simu_fit)$v0[1,]
which(abs(v0s) < 0.01)
which(v0s > 1.2)
which(v0s < -1)
# Neutral 121, Positive 197, Negative 24
