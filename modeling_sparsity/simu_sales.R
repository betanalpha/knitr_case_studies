############################################################
# Initial setup
############################################################

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

util <- new.env()
source('stan_utility.R', local=util)
source('sparse_utility.R', local=util)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_trans <- c("#DCBCBC80")
c_dark_trans <- c("#8F272780")
c_green_trans <- c("#00FF0080")

c_light_teal <- c("#6B8E8E")
c_mid_teal <- c("#487575")
c_dark_teal <- c("#1D4F4F")

par(family="CMU Serif", las=1, bty="l", cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 5))

############################################################
# Create data
############################################################

set.seed(593838393)

N_obs_days <- 370

obs_day_idx <- unlist(lapply(1:N_obs_days, function(d) rep(d, sample(c(3, 4, 5), 1))))

N_obs <- length(obs_day_idx)

N_pred_days <- 100
pred_day_idx <- 1:N_pred_days + N_obs_days + 30

simu <- stan(file='stan_programs/simu_sales.stan', 
             data = list("N_obs" = N_obs, 
                         "N_obs_days" = N_obs_days, "obs_day_idx" = obs_day_idx,
                         "N_pred_days" = N_pred_days, "pred_day_idx" = pred_day_idx), 
             iter=1, chains=1, seed=2948399, algorithm="Fixed_param")

y_obs <- extract(simu)$y_obs[1,]
y_pred <- extract(simu)$y_pred[1,]

stan_rdump(c("y_obs", "N_obs", "N_obs_days", "obs_day_idx",
             "y_pred", "N_pred_days", "pred_day_idx"),
           file="data/sales.data.R")

alpha <- log(226)
delta_season <- log(1.3)
phi <- 10
inner_tau_day <- log(1.005)
outer_tau_day <- log(1.25)
gamma <- 0.9
inv_psi <-  0.001

delta_days <- extract(simu)$delta_days[1,]

stan_rdump(c("alpha", "delta_season", "phi", 
             "inner_tau_day", "outer_tau_day", 
             "gamma", "inv_psi", "delta_days"),
           file="data/sales.truth.R")