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

############################################################
# Explore data
############################################################

data <- read_rdump('data/sales.data.R')

plot(c(data$obs_day_idx, data$pred_day_idx),
     c(data$y_obs, data$y_pred),
     col="black", pch=16,
     xlab="Day", xlim=c(0, 500),
     ylab="Total Sales", ylim=c(0, 400))

truth <- read_rdump('data/sales.truth.R')

############################################################
# Normal Population Model
############################################################

# Non-center all of the day excess parameters
data$cp_idx <- vector()
data$K_cp <- length(data$cp_idx)

data$ncp_idx <- setdiff(1:data$N_obs_days, data$cp_idx)
data$K_ncp <- length(data$ncp_idx)

normal_fit <- stan(file='stan_programs/fit_sales_normal.stan', data=data, 
                   seed=4938483, refresh=1000)

# Check diagnostics
util$check_all_diagnostics(normal_fit)

# Check marginal posteriors
samples = extract(normal_fit)

par(mfrow=c(2, 3))

hist(samples$alpha, breaks=seq(5.40, 5.45, 0.05 / 50),
     main="", xlab="alpha", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$alpha, col=c_light, lty=1, lw=2)

hist(samples$delta_season, seq(0.2, 0.3, 0.1 / 50),
     main="", xlab="delta_season", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$delta_season, col=c_light, lty=1, lw=2)

hist(samples$phi, seq(9.8, 10.2, 0.4 / 50), 
     main="", xlab="phi", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$phi, col=c_light, lty=1, lw=2)

hist(samples$tau_day, seq(0, 0.1, 0.1 / 50),
     main="", xlab="tau_day", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(samples$milli_inv_psi, seq(0, 3, 3 / 50),
     main="", xlab="milli_inv_psi", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=1000 * truth$inv_psi, col=c_light, lty=1, lw=2)

par(mfrow=c(2, 1))

util$plot_excess_variation(samples, data, truth, "Normal Population Model")

util$plot_excess_variation_residual(samples, data, truth, "Normal Population Model")

infer_crps_comp <- data.frame("Normal" = util$sum_empirical_crps(samples$delta_days, truth$delta_days))
row.names(infer_crps_comp) = "Inferential ECRPS"

print(infer_crps_comp)

# Check expected sales time series
par(mfrow=c(2, 1))

util$plot_expected_sales(samples, data, "Normal Population Model")

util$plot_expected_sales_residual(samples, data, "Normal Population Model")

# Check posterior retrodictions and predictions
util$plot_dictions(samples, data, "Normal Population Model")

util$plot_dictions_residual(samples, data, "Normal Population Model")

# Evaluate predictive performance with cumulative risk probability score
pred_crps_comp <- data.frame("Normal" = util$sum_empirical_crps(samples$y_post_pre, data$y_pred))
row.names(pred_crps_comp) = "Predictive ECRPS"

print(pred_crps_comp)

############################################################
# Laplace Population Model
############################################################

data$cp_idx <- which(abs(sapply(1:data$N_obs_day, 
                                function(n) mean(samples$delta_days[, n]))) > 0.01)
data$K_cp <- length(data$cp_idx)

data$ncp_idx <- setdiff(1:data$N_obs_days, data$cp_idx)
data$K_ncp <- length(data$ncp_idx)

laplace_fit <- stan(file='stan_programs/fit_sales_laplace.stan', data=data, 
                    seed=4938483, refresh=1000)

# Check diagnostics
util$check_all_diagnostics(laplace_fit)

# Check marginal posteriors
samples = extract(laplace_fit)

par(mfrow=c(2, 3))

hist(samples$alpha, breaks=seq(5.40, 5.45, 0.05 / 50),
     main="", xlab="alpha", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$alpha, col=c_light, lty=1, lw=2)

hist(samples$delta_season, seq(0.2, 0.3, 0.1 / 50),
     main="", xlab="delta_season", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$delta_season, col=c_light, lty=1, lw=2)

hist(samples$phi, seq(9.8, 10.2, 0.4 / 50), 
     main="", xlab="phi", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$phi, col=c_light, lty=1, lw=2)

hist(samples$tau_day, seq(0, 0.1, 0.1 / 50),
     main="", xlab="tau_day", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(samples$milli_inv_psi, seq(0, 3, 3 / 50),
     main="", xlab="milli_inv_psi", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=1000 * truth$inv_psi, col=c_light, lty=1, lw=2)

par(mfrow=c(2, 1))

util$plot_excess_variation(samples, data, truth, "Laplace Population Model")

util$plot_excess_variation_residual(samples, data, truth, "Laplace Population Model")

infer_crps_comp["Laplace"] = util$sum_empirical_crps(samples$delta_days, truth$delta_days)
print(infer_crps_comp)

# Check expected sales time series
par(mfrow=c(2, 1))

util$plot_expected_sales(samples, data, "Laplace Population Model")

util$plot_expected_sales_residual(samples, data, "Laplace Population Model")

# Check posterior retrodictions and predictions
util$plot_dictions(samples, data, "Laplace Population Model")

util$plot_dictions_residual(samples, data, "Laplace Population Model")

# Evaluate predictive performance with cumulative risk probability score
pred_crps_comp["Laplace"] = util$sum_empirical_rps(samples$y_post_pre, data$y_pred)
print(pred_crps_comp)

############################################################
# Cauchy Population Model
############################################################

set.seed(94858292)

samples <- extract(normal_fit)

data$cp_idx <- which(abs(sapply(1:data$N_obs_day, 
                                function(n) mean(samples$delta_days[, n]))) > 0.05)
data$K_cp <- length(data$cp_idx)

data$ncp_idx <- setdiff(1:data$N_obs_days, data$cp_idx)
data$K_ncp <- length(data$ncp_idx)

# Using same non-centering as before

cauchy_fit <- stan(file='stan_programs/fit_sales_cauchy.stan', data=data, 
                    seed=4938483, refresh=1000, init_r=1)

# Check diagnostics
util$check_all_diagnostics(cauchy_fit)

sapply(1:4, function(c) get_sampler_params(cauchy_fit, inc_warmup=FALSE)[[c]][,'stepsize__'][1])

# Check marginal posteriors
samples = extract(cauchy_fit)

par(mfrow=c(2, 3))

hist(samples$alpha, breaks=seq(5.40, 5.45, 0.05 / 50),
     main="", xlab="alpha", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$alpha, col=c_light, lty=1, lw=2)

hist(samples$delta_season, seq(0.2, 0.3, 0.1 / 50),
     main="", xlab="delta_season", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$delta_season, col=c_light, lty=1, lw=2)

hist(samples$phi, seq(0, 15, 15 / 50), 
     main="", xlab="phi", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$phi, col=c_light, lty=1, lw=2)

hist(samples$tau_day, seq(0, 0.1, 0.1 / 50),
     main="", xlab="tau_day", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(samples$milli_inv_psi, seq(0, 50, 50 / 50),
     main="", xlab="milli_inv_psi", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=1000 * truth$inv_psi, col=c_light, lty=1, lw=2)

par(mfrow=c(2, 1))

util$plot_excess_variation(samples, data, truth, "Cauchy Population Model")

util$plot_excess_variation_residual(samples, data, truth, "Cauchy Population Model")

infer_crps_comp["Cauchy"] = util$sum_empirical_crps(samples$delta_days, truth$delta_days)
print(infer_crps_comp)

# Check expected sales time series
par(mfrow=c(2, 1))

util$plot_expected_sales(samples, data, "Cauchy Population Model")

util$plot_expected_sales_residual(samples, data, "Cauchy Population Model")

# Check posterior retrodictions and predictions
util$plot_dictions(samples, data, "Cauchy Population Model")

util$plot_dictions_residual(samples, data, "Cauchy Population Model")

# Evaluate predictive performance with cumulative risk probability score
pred_crps_comp["Cauchy"] = util$sum_empirical_rps(samples$y_post_pre, data$y_pred)
print(pred_crps_comp)

############################################################
# Binary Normal Mixture Population Model
############################################################

data$w <- rep(0.25, N_obs_days)
data$w[data$cp_idx] <- 1

data$w <- c(0.25, 0.25, 0.25, 0.25, 1, 0.25, 1, 0.25, 0.25, 0.25, 1, 0.25, 
  0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 
  1, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
  0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
  0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 
  0.25, 0.25, 1, 0.25, 1, 0.25, 0.25, 0.25, 1, 1, 0.25, 1, 0.25, 
  0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
  0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 
  0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 
  0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
  0.25, 0.25, 1, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 
  0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
  0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 0.25, 0.25, 
  0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 
  0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 
  0.25, 1, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
  0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
  1, 0.25, 1, 1, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 
  0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
  0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
  0.25, 1, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 0.25, 0.25, 
  0.25, 1, 0.25, 0.25, 1, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 0.25, 
  0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
  0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
  0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 0.25, 0.25, 
  1, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 
  1, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 1, 1, 0.25, 0.25, 
  0.25, 0.25, 0.25, 0.25, 1, 1, 0.25, 1, 1, 0.25, 0.25, 0.25, 0.25, 
  0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 0.25, 0.25, 
  0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 
  0.25, 0.25, 1, 0.25, 0.25, 0.25, 1, 0.25, 1, 0.25, 0.25, 0.25, 
  0.25, 1, 0.25, 0.25, 1, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
  0.25, 0.25)

.Random.seed <- c(10403L, 40L, -809302514L, -1649027925L, 836617953L, -1249118417L, 
                  1278836891L, 1685176797L, -1654588936L, 1905218982L, 1418063728L, 
                  -900832307L, 10935867L, 223600251L, 918001093L, -193246870L, 
                  -715013192L, -1147444885L, -48170971L, 287154586L, 898872228L, 
                  1985097474L, 813880538L, -1462121185L, 464361994L, -1679288726L, 
                  -1993959960L, 1986444395L, -1645791546L, 1218736261L, 104697543L, 
                  487323447L, 524642538L, -1186938912L, -1459728208L, 1573319732L, 
                  -922967540L, 496703551L, -2005034741L, -1162069396L, -588094218L, 
                  2084118371L, 1336334429L, 42159845L, -367497111L, -758974312L, 
                  -123439497L, -1445243774L, -1759287408L, -2090337883L, -1207465398L, 
                  512891L, -778291492L, 1508498582L, -1911834872L, 335919586L, 
                  -601280334L, -519957804L, 206864330L, -403042822L, 1769430502L, 
                  966491242L, -2142052785L, 639379867L, 964716993L, 1450175052L, 
                  -365018913L, -1209520458L, 1126038731L, -1197886125L, 1722141607L, 
                  -1199047730L, 1681976421L, -1876482649L, -425952786L, -834685866L, 
                  1730983462L, -1374057795L, -1300668353L, 2095553478L, -409587153L, 
                  1041759712L, -1505991315L, -1809682093L, 393386179L, 286872881L, 
                  1062790857L, -482369449L, -1490862928L, 1599068955L, 2077707579L, 
                  -218063015L, 556520195L, 117761664L, -815357251L, -153321107L, 
                  1395457942L, -418515380L, 1741531716L, 879268998L, 618519779L, 
                  665619798L, 952670250L, 297368708L, 1441631573L, -312543908L, 
                  -444007866L, -2084615288L, -1594160057L, -208111664L, 943488324L, 
                  1643677429L, -1442519946L, 79931351L, 1789585057L, 1008119529L, 
                  1477522689L, 959184526L, 1484627404L, -655905929L, 766565263L, 
                  -1800491409L, -1436807012L, 1902683776L, -1161109379L, 376785194L, 
                  1865022879L, 1645315816L, 1719529771L, 1302645265L, 2102266452L, 
                  -1965434073L, -120290649L, -118956406L, 1267690685L, -1001026116L, 
                  159582074L, 911065087L, 1679083262L, 1567882608L, 1848371673L, 
                  -1016171307L, 718969749L, 1134295365L, 271446877L, 2009938777L, 
                  -92941395L, -556135669L, -1849897848L, 1617656258L, 2053076467L, 
                  -1435652697L, 681414953L, 1993110684L, 1970802862L, -902581651L, 
                  -914370439L, -2122194650L, 1181431448L, -862623383L, -500929858L, 
                  -768624268L, 614616720L, 262822418L, 55117829L, 793387909L, -471556996L, 
                  -1617674827L, 926192080L, 2131335012L, -13819516L, -103189180L, 
                  -421002356L, 1168301701L, -1751731071L, -1573154676L, 1813326885L, 
                  -387334501L, 1112159171L, -220128778L, 1975424258L, 1264689871L, 
                  -997124153L, -2089682654L, -1623245438L, -1290180552L, -1221412066L, 
                  -778777468L, 2138238340L, -646798295L, -1239262575L, 552796697L, 
                  1348616574L, -144987778L, 1258715279L, -1936520233L, 2045187265L, 
                  2049723351L, -162145452L, 1479880025L, 1384696419L, 191385880L, 
                  1865399051L, -956415030L, 1582508754L, 1918936458L, -768794464L, 
                  -96937557L, 1396363941L, -394113433L, 765288709L, 1583549289L, 
                  2093964499L, 1643682059L, 2030422695L, -561671778L, 1806106102L, 
                  -636217839L, 778744373L, -1303863103L, 159629913L, 815279723L, 
                  -1324586728L, 1022926588L, -1999268082L, 1181263850L, -649907656L, 
                  -1256689846L, 1464015558L, -1184062827L, -223744786L, 2106675482L, 
                  858663020L, -1032182637L, -1808803529L, 1311778326L, -1139278779L, 
                  1721839891L, 1157423655L, 1000192085L, -1836182789L, -206233435L, 
                  -180870101L, 1716309376L, 989031072L, 1790442471L, 639569958L, 
                  384038419L, -977917626L, -792580150L, 1841803693L, -173293031L, 
                  -411274345L, -1776500492L, 1377798778L, -1077106283L, -349209973L, 
                  -774490401L, -1233572039L, -683888514L, -59702360L, 1407445834L, 
                  -949296052L, -980409565L, -1734407977L, 1141361783L, -41617754L, 
                  -1628720697L, 1403626520L, -889591515L, -262662905L, -1143335352L, 
                  116111341L, 701702832L, 2047137284L, 623236480L, 1093741154L, 
                  392118927L, -74217547L, 1025636891L, -142625745L, -1099822570L, 
                  1548730617L, 1965207878L, -216175698L, -728665762L, 517539541L, 
                  35751103L, 243980891L, -705264552L, -939374650L, -275884252L, 
                  939871047L, -495109493L, -1906814392L, 1323393329L, -105022075L, 
                  -1365165508L, -2056433065L, 467307105L, -2051658777L, 470078325L, 
                  -1133749909L, 2040299753L, 1177484256L, 1841931556L, -1676819337L, 
                  -911451489L, -1157751754L, 835177681L, 1643897980L, 1527957938L, 
                  2348978L, -881387386L, -1507999734L, -612486017L, 46847612L, 
                  -1799422625L, 2024118669L, -312907091L, -726467362L, 203820618L, 
                  1026147730L, 1989244774L, 1671372849L, 1878891750L, 1196138206L, 
                  622469620L, -1267854396L, 1392467275L, -56696252L, -853923739L, 
                  -2000398591L, -740453798L, 1246112720L, -794673496L, 1442680937L, 
                  783655076L, 700013672L, -1505782159L, 994550737L, -418180376L, 
                  -1966920885L, -365778063L, 488365994L, -556628251L, -2127004524L, 
                  -386401573L, 1745451968L, 1205113358L, 432773061L, -2117839630L, 
                  -557726508L, 1482643433L, -1651156357L, 1847675515L, 1316770L, 
                  1790574982L, 1854834157L, -1473940789L, 1794678962L, -355817063L, 
                  1217133699L, 1891774341L, -17942259L, 61987390L, -2129214625L, 
                  -846810315L, 376893477L, -1073314196L, 1262161900L, 717193001L, 
                  1100479742L, -1829011200L, -139601754L, -77521541L, -501967070L, 
                  1257155013L, -802734293L, -2098270454L, 827861875L, 1726788307L, 
                  -281554663L, 49888470L, 637741180L, -499082275L, -1217164996L, 
                  668990235L, -815615977L, 1945850457L, 1741817277L, 1851548601L, 
                  870124805L, -1439787676L, 1166421734L, -868138993L, 921203235L, 
                  1210522568L, -75112018L, 1304013213L, -872295368L, 1077866087L, 
                  1682260573L, -2108095741L, 520235509L, 1162003169L, -225664611L, 
                  -1372071407L, -2112165590L, -1979484341L, -1630877256L, -636848109L, 
                  1431591212L, -2123795028L, 843974784L, -532877801L, 953393990L, 
                  819740317L, 708479931L, 467554348L, -1227368880L, -999625313L, 
                  1605458723L, -881784612L, -1157183709L, 87900280L, 1573180896L, 
                  -1998122005L, 1755844923L, -293806485L, -1098544874L, -137597504L, 
                  -1222868541L, -590059898L, -2098292261L, -1524293286L, -120096091L, 
                  651930973L, 1550266628L, -29876983L, -1914481919L, -1314782060L, 
                  -2079006477L, 1104688370L, -1113441369L, -1628981824L, -1031558852L, 
                  424874107L, -1292723418L, 2135840786L, 58420209L, 534241037L, 
                  -841478802L, 429792318L, 203584717L, 1176050954L, 1669982184L, 
                  -1601595454L, 285923460L, 1164713879L, -699533839L, -1814358291L, 
                  -904191595L, -1649752477L, 1042742332L, 498534036L, -478231205L, 
                  -873795188L, -1427797050L, 1402337481L, 1702891913L, -1731273055L, 
                  762039641L, -933009058L, -889006415L, -1487199978L, -1769608932L, 
                  -1964574055L, -1457647348L, -82195063L, 342047643L, 127328441L, 
                  -1111107843L, -1316575983L, 1935500293L, 693429999L, 813124143L, 
                  -1470564320L, 1688073866L, 1894566102L, -759839983L, 892945585L, 
                  1743117242L, 226105754L, -1131174359L, -1641596741L, 1616028097L, 
                  143536006L, -360106430L, -219779721L, 1878453250L, 881461792L, 
                  1792173365L, -1460229752L, -1256843033L, 181839880L, -1714444504L, 
                  -908170874L, -652123112L, -221576806L, 1893421092L, -1337017914L, 
                  1810472933L, 1129781333L, 1414512861L, -1383904912L, 1553945052L, 
                  1959714339L, -737463321L, -1208684064L, -1806682705L, 274798480L, 
                  283185078L, 1394174670L, 1056174018L, 476966481L, 992957759L, 
                  1314393149L, 1626051874L, -1959865212L, 1934035236L, -1578656572L, 
                  1388860634L, -385346117L, -797213886L, -265885879L, -1369523083L, 
                  -109535459L, 63446241L, 142901344L, -869233916L, 1649283766L, 
                  -156001219L, -972161338L, -1785475312L, -2111545868L, -672239067L, 
                  -1793021768L, -1829083699L, 212149810L, -408582437L, 4398326L, 
                  -815423659L, 928236454L, 1558293485L, -1294777372L, 1687563396L, 
                  -1112158858L, 401013199L, 1118621393L, 1457662597L, 1758915456L, 
                  -1972203477L, 1616310339L, 568364066L, -2142105725L, 1113152765L, 
                  -880727309L, 463725441L, 68464333L, -1197337398L, 1535690804L, 
                  -336351145L, -1449160593L, -1777629838L, -193380907L, 128629040L, 
                  327722364L, 596765589L, -943329089L, 1178566603L, -1871500841L, 
                  -1693343807L, -648319107L, -4347453L, 258884252L, -1520999381L, 
                  -1328082201L, 1783130926L, 1423617136L, 940552L, 1854841677L, 
                  -808321327L, 1280062860L, 1712918026L, -1013656639L, 807777202L, 
                  -1454957880L, -1823314840L, 1701357697L, 539651538L, 193736663L, 
                  1146364944L, -1291420497L, -1651535721L, -1799524634L, 854016906L, 
                  -53837673L, -125343026L, 1174705402L, -621750140L, -953554677L, 
                  -1079572619L, 497705036L, -1834145204L, 1667496812L, 972695281L, 
                  230385907L, -496694749L, -1199420939L, -11293146L, 323874436L, 
                  -1462133397L, 879957117L, 2095594610L, -1671142743L)

data$w == 0.25

normal_mixture_fit <- stan(file='stan_programs/fit_sales_normal_mixture.stan',
                           data=data, seed=4938483, refresh=1000)

util$check_all_diagnostics(normal_mixture_fit)

partition <- util$partition_div(normal_mixture_fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

par(mfrow=c(1, 3))

for (k in c(36, 238, 311)) {
  name_x <- paste("delta_days_tilde[", k, "]", sep='')
  name_y <- "inner_tau_day"
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("Day ", k),
       xlab=name_x, xlim=c(-10, 10), 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-7, -3))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
}

unpermuted_samples <- extract(normal_mixture_fit, permute=FALSE)

par(mfrow=c(2, 2))

for (k in c(36, 238, 311))
  for (c in 1:4)
    plot(1:1000, unpermuted_samples[, c, k], type="l", lwd=1, col=c_dark,
         main=paste("Chain", c, sep=""),
         xlab="Iteration",  xlim=c(1, 1000),
         ylab=paste("delta_days_tilde[", k, "]", sep=""), ylim=c(-10, 10))






data$cp_idx <- which(abs(sapply(1:data$N_obs_day, 
                                function(n) mean(extract(normal_fit)$delta_days[, n]))) > 0.025)
data$K_cp <- length(data$cp_idx)

data$ncp_idx <- setdiff(1:data$N_obs_days, data$cp_idx)
data$K_ncp <- length(data$ncp_idx)


data$cp_idx <- which(abs(sapply(1:data$N_obs_day, 
                                function(n) mean(extract(normal_fit)$delta_days[, n]))) > 0.05)
data$K_cp <- length(data$cp_idx)

data$ncp_idx <- setdiff(1:data$N_obs_days, data$cp_idx)
data$K_ncp <- length(data$ncp_idx)

# Note the same (non)centering pattern is still being used
set.seed(83855834)

normal_mixture_fit <- stan(file='stan_programs/fit_sales_normal_mixture.stan',
                           data=data, seed=4938483, refresh=1000, init_r=0.1)

# Check diagnostics
util$check_all_diagnostics(normal_mixture_fit)

sapply(1:4, function(c) get_sampler_params(normal_mixture_fit, inc_warmup=FALSE)[[c]][,'stepsize__'][1])

unpermuted_samples <- extract(normal_mixture_fit, permute=FALSE)

par(mfrow=c(2, 2))

plot_idx <- data$K_cp + 20

for (c in 1:4) {
  plot(1:1000, unpermuted_samples[,c,plot_idx], type="l", lwd=1, col=c_dark,
       main=paste("Chain", c),
       xlab="Iteration",  xlim=c(1, 1000),
       ylab="delta_days_ncp[20]", ylim=c(-40, 40))
}

par(mfrow=c(2, 2))

day_idx <- data$ncp_idx[20]
plot_idx <- data$N_obs_days + 8 + day_idx

for (c in 1:4) {
  plot(1:1000, unpermuted_samples[,c,plot_idx], type="l", lwd=1, col=c_dark,
       main=paste("Chain", c),
       xlab="Iteration",  xlim=c(1, 1000),
       ylab=paste("delta_days[", day_idx, "]", sep=""), ylim=c(-0.4, 0.4))
  
  abline(h=-truth$inner_tau_day, col=c_light, lty=1, lw=1)
  abline(h=+truth$inner_tau_day, col=c_light, lty=1, lw=1)
  
  abline(h=-truth$outer_tau_day, col=c_mid, lty=1, lw=1)
  abline(h=+truth$outer_tau_day, col=c_mid, lty=1, lw=1)
}

# Truth value looks okay
truth$delta_days[day_idx]

# Expected sales looks okay
mu <- exp(truth$alpha + truth$delta_season * sin(2 * pi * 25 / 360 + truth$phi) + truth$delta_days[25])
print(mu)

# What about the data?
data$y_obs[which(data$obs_day_idx == day_idx)]

# Ah, the data all fluctuated low, so the likelihood function for this day
# concentrates at smaller values that will influence the corresponding 
# delta_day away from zero, but not quite far enough away for the posterior
# distribution to decouple from the inner core.

# The relatively sharp transition between the inner and outer components 
# can be difficult to navigate, resulting in this weak metastable behavior 
# and slow exploration.

# We can smooth out this interface by using a slightly more heavy-tailed
# model for the inner core.


mixture_fit <- stan(file='stan_programs/fit_sales_student_normal_mixture.stan',
                           data=data, seed=74938483, refresh=1000)

# Check diagnostics
util$check_all_diagnostics(mixture_fit)

unpermuted_samples <- extract(mixture_fit, permute=FALSE)

par(mfrow=c(2, 2))

plot_idx <- data$K_cp + 2

for (c in 1:4) {
  plot(1:1000, unpermuted_samples[,c,plot_idx], type="l", lwd=1, col=c_dark,
       main=paste("Chain", c),
       xlab="Iteration",  xlim=c(1, 1000),
       ylab="delta_days_ncp[20]", ylim=c(-40, 40))
}

data$ncp_idx[2]
data$ncp_idx[36]

par(mfrow=c(2, 2))

day_idx <- data$ncp_idx[2]
plot_idx <- data$N_obs_days + 8 + day_idx

for (c in 1:4) {
  plot(1:1000, unpermuted_samples[,c,plot_idx], type="l", lwd=1, col=c_dark,
       main=paste("Chain", c),
       xlab="Iteration",  xlim=c(1, 1000),
       ylab=paste("delta_days[", day_idx, "]", sep=""), ylim=c(-0.4, 0.4))
  
  abline(h=-truth$inner_tau_day, col=c_light, lty=1, lw=1)
  abline(h=+truth$inner_tau_day, col=c_light, lty=1, lw=1)
  
  abline(h=-truth$outer_tau_day, col=c_mid, lty=1, lw=1)
  abline(h=+truth$outer_tau_day, col=c_mid, lty=1, lw=1)
}

# Truth value looks okay
truth$delta_days[day_idx]

# Expected sales looks okay
mu <- exp(truth$alpha + truth$delta_season * sin(2 * pi * 25 / 360 + truth$phi) + truth$delta_days[25])
print(mu)

# What about the data?
data$y_obs[which(data$obs_day_idx == day_idx)]

# Check marginal posteriors
samples = extract(mixture_fit)

par(mfrow=c(2, 4))

hist(samples$alpha, breaks=seq(5.40, 5.45, 0.05 / 50),
     main="", xlab="alpha", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$alpha, col=c_light, lty=1, lw=2)

hist(samples$delta_season, seq(0.2, 0.3, 0.1 / 50),
     main="", xlab="delta_season", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$delta_season, col=c_light, lty=1, lw=2)

hist(samples$phi, seq(9.8, 10.2, 0.4 / 50), 
     main="", xlab="phi", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$phi, col=c_light, lty=1, lw=2)

hist(samples$inner_tau_day, seq(0, 0.1, 0.1 / 50),
     main="", xlab="inner_tau_day", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$inner_tau_day, col=c_light, lty=1, lw=2)

hist(samples$outer_tau_day, seq(0, 0.5, 0.5 / 50),
     main="", xlab="outer_tau_day", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$outer_tau_day, col=c_light, lty=1, lw=2)

hist(samples$gamma, seq(0.75, 1.0, 0.25 / 50),
     main="", xlab="gamma", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$gamma, col=c_light, lty=1, lw=2)

hist(samples$milli_inv_psi, seq(0, 3, 3 / 50),
     main="", xlab="milli_inv_psi", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=1000 * truth$inv_psi, col=c_light, lty=1, lw=2)

par(mfrow=c(2, 1))

util$plot_excess_variation(samples, data, truth, "Normal Mixture Population Model")

util$plot_excess_variation_residual(samples, data, truth, "Normal Mixture Population Model")

infer_crps_comp["Normal Mixture"] = util$sum_empirical_crps(samples$delta_days, truth$delta_days)
print(infer_crps_comp)

# Check expected sales time series
par(mfrow=c(2, 1))

util$plot_expected_sales(samples, data, "Normal Mixture Population Model")

util$plot_expected_sales_residual(samples, data, "Normal Mixture Population Model")

# Check posterior retrodictions and predictions
util$plot_dictions(samples, data, "Normal Mixture Population Model")

util$plot_dictions_residual(samples, data, "Normal Mixture Population Model")

# Evaluate predictive performance with cumulative risk probability score
pred_crps_comp["Normal Mixture"] = util$sum_empirical_rps(samples$y_post_pre, data$y_pred)
print(pred_crps_comp)

############################################################
# Horseshoe Population Model
############################################################

test_seed <- .Random.seed

rnorm(5, 0, 1)
rnorm(5, 0, 1)

.Random.seed <- test_seed

rnorm(5, 0, 1)

dput(.Random.seed)


.Random.seed <- c(10403L, 37L, -809302514L, -1649027925L, 836617953L, -1249118417L, 
                  1278836891L, 1685176797L, -1654588936L, 1905218982L, 1418063728L, 
                  -900832307L, 10935867L, 223600251L, 918001093L, -193246870L, 
                  -715013192L, -1147444885L, -48170971L, 287154586L, 898872228L, 
                  1985097474L, 813880538L, -1462121185L, 464361994L, -1679288726L, 
                  -1993959960L, 1986444395L, -1645791546L, 1218736261L, 104697543L, 
                  487323447L, 524642538L, -1186938912L, -1459728208L, 1573319732L, 
                  -922967540L, 496703551L, -2005034741L, -1162069396L, -588094218L, 
                  2084118371L, 1336334429L, 42159845L, -367497111L, -758974312L, 
                  -123439497L, -1445243774L, -1759287408L, -2090337883L, -1207465398L, 
                  512891L, -778291492L, 1508498582L, -1911834872L, 335919586L, 
                  -601280334L, -519957804L, 206864330L, -403042822L, 1769430502L, 
                  966491242L, -2142052785L, 639379867L, 964716993L, 1450175052L, 
                  -365018913L, -1209520458L, 1126038731L, -1197886125L, 1722141607L, 
                  -1199047730L, 1681976421L, -1876482649L, -425952786L, -834685866L, 
                  1730983462L, -1374057795L, -1300668353L, 2095553478L, -409587153L, 
                  1041759712L, -1505991315L, -1809682093L, 393386179L, 286872881L, 
                  1062790857L, -482369449L, -1490862928L, 1599068955L, 2077707579L, 
                  -218063015L, 556520195L, 117761664L, -815357251L, -153321107L, 
                  1395457942L, -418515380L, 1741531716L, 879268998L, 618519779L, 
                  665619798L, 952670250L, 297368708L, 1441631573L, -312543908L, 
                  -444007866L, -2084615288L, -1594160057L, -208111664L, 943488324L, 
                  1643677429L, -1442519946L, 79931351L, 1789585057L, 1008119529L, 
                  1477522689L, 959184526L, 1484627404L, -655905929L, 766565263L, 
                  -1800491409L, -1436807012L, 1902683776L, -1161109379L, 376785194L, 
                  1865022879L, 1645315816L, 1719529771L, 1302645265L, 2102266452L, 
                  -1965434073L, -120290649L, -118956406L, 1267690685L, -1001026116L, 
                  159582074L, 911065087L, 1679083262L, 1567882608L, 1848371673L, 
                  -1016171307L, 718969749L, 1134295365L, 271446877L, 2009938777L, 
                  -92941395L, -556135669L, -1849897848L, 1617656258L, 2053076467L, 
                  -1435652697L, 681414953L, 1993110684L, 1970802862L, -902581651L, 
                  -914370439L, -2122194650L, 1181431448L, -862623383L, -500929858L, 
                  -768624268L, 614616720L, 262822418L, 55117829L, 793387909L, -471556996L, 
                  -1617674827L, 926192080L, 2131335012L, -13819516L, -103189180L, 
                  -421002356L, 1168301701L, -1751731071L, -1573154676L, 1813326885L, 
                  -387334501L, 1112159171L, -220128778L, 1975424258L, 1264689871L, 
                  -997124153L, -2089682654L, -1623245438L, -1290180552L, -1221412066L, 
                  -778777468L, 2138238340L, -646798295L, -1239262575L, 552796697L, 
                  1348616574L, -144987778L, 1258715279L, -1936520233L, 2045187265L, 
                  2049723351L, -162145452L, 1479880025L, 1384696419L, 191385880L, 
                  1865399051L, -956415030L, 1582508754L, 1918936458L, -768794464L, 
                  -96937557L, 1396363941L, -394113433L, 765288709L, 1583549289L, 
                  2093964499L, 1643682059L, 2030422695L, -561671778L, 1806106102L, 
                  -636217839L, 778744373L, -1303863103L, 159629913L, 815279723L, 
                  -1324586728L, 1022926588L, -1999268082L, 1181263850L, -649907656L, 
                  -1256689846L, 1464015558L, -1184062827L, -223744786L, 2106675482L, 
                  858663020L, -1032182637L, -1808803529L, 1311778326L, -1139278779L, 
                  1721839891L, 1157423655L, 1000192085L, -1836182789L, -206233435L, 
                  -180870101L, 1716309376L, 989031072L, 1790442471L, 639569958L, 
                  384038419L, -977917626L, -792580150L, 1841803693L, -173293031L, 
                  -411274345L, -1776500492L, 1377798778L, -1077106283L, -349209973L, 
                  -774490401L, -1233572039L, -683888514L, -59702360L, 1407445834L, 
                  -949296052L, -980409565L, -1734407977L, 1141361783L, -41617754L, 
                  -1628720697L, 1403626520L, -889591515L, -262662905L, -1143335352L, 
                  116111341L, 701702832L, 2047137284L, 623236480L, 1093741154L, 
                  392118927L, -74217547L, 1025636891L, -142625745L, -1099822570L, 
                  1548730617L, 1965207878L, -216175698L, -728665762L, 517539541L, 
                  35751103L, 243980891L, -705264552L, -939374650L, -275884252L, 
                  939871047L, -495109493L, -1906814392L, 1323393329L, -105022075L, 
                  -1365165508L, -2056433065L, 467307105L, -2051658777L, 470078325L, 
                  -1133749909L, 2040299753L, 1177484256L, 1841931556L, -1676819337L, 
                  -911451489L, -1157751754L, 835177681L, 1643897980L, 1527957938L, 
                  2348978L, -881387386L, -1507999734L, -612486017L, 46847612L, 
                  -1799422625L, 2024118669L, -312907091L, -726467362L, 203820618L, 
                  1026147730L, 1989244774L, 1671372849L, 1878891750L, 1196138206L, 
                  622469620L, -1267854396L, 1392467275L, -56696252L, -853923739L, 
                  -2000398591L, -740453798L, 1246112720L, -794673496L, 1442680937L, 
                  783655076L, 700013672L, -1505782159L, 994550737L, -418180376L, 
                  -1966920885L, -365778063L, 488365994L, -556628251L, -2127004524L, 
                  -386401573L, 1745451968L, 1205113358L, 432773061L, -2117839630L, 
                  -557726508L, 1482643433L, -1651156357L, 1847675515L, 1316770L, 
                  1790574982L, 1854834157L, -1473940789L, 1794678962L, -355817063L, 
                  1217133699L, 1891774341L, -17942259L, 61987390L, -2129214625L, 
                  -846810315L, 376893477L, -1073314196L, 1262161900L, 717193001L, 
                  1100479742L, -1829011200L, -139601754L, -77521541L, -501967070L, 
                  1257155013L, -802734293L, -2098270454L, 827861875L, 1726788307L, 
                  -281554663L, 49888470L, 637741180L, -499082275L, -1217164996L, 
                  668990235L, -815615977L, 1945850457L, 1741817277L, 1851548601L, 
                  870124805L, -1439787676L, 1166421734L, -868138993L, 921203235L, 
                  1210522568L, -75112018L, 1304013213L, -872295368L, 1077866087L, 
                  1682260573L, -2108095741L, 520235509L, 1162003169L, -225664611L, 
                  -1372071407L, -2112165590L, -1979484341L, -1630877256L, -636848109L, 
                  1431591212L, -2123795028L, 843974784L, -532877801L, 953393990L, 
                  819740317L, 708479931L, 467554348L, -1227368880L, -999625313L, 
                  1605458723L, -881784612L, -1157183709L, 87900280L, 1573180896L, 
                  -1998122005L, 1755844923L, -293806485L, -1098544874L, -137597504L, 
                  -1222868541L, -590059898L, -2098292261L, -1524293286L, -120096091L, 
                  651930973L, 1550266628L, -29876983L, -1914481919L, -1314782060L, 
                  -2079006477L, 1104688370L, -1113441369L, -1628981824L, -1031558852L, 
                  424874107L, -1292723418L, 2135840786L, 58420209L, 534241037L, 
                  -841478802L, 429792318L, 203584717L, 1176050954L, 1669982184L, 
                  -1601595454L, 285923460L, 1164713879L, -699533839L, -1814358291L, 
                  -904191595L, -1649752477L, 1042742332L, 498534036L, -478231205L, 
                  -873795188L, -1427797050L, 1402337481L, 1702891913L, -1731273055L, 
                  762039641L, -933009058L, -889006415L, -1487199978L, -1769608932L, 
                  -1964574055L, -1457647348L, -82195063L, 342047643L, 127328441L, 
                  -1111107843L, -1316575983L, 1935500293L, 693429999L, 813124143L, 
                  -1470564320L, 1688073866L, 1894566102L, -759839983L, 892945585L, 
                  1743117242L, 226105754L, -1131174359L, -1641596741L, 1616028097L, 
                  143536006L, -360106430L, -219779721L, 1878453250L, 881461792L, 
                  1792173365L, -1460229752L, -1256843033L, 181839880L, -1714444504L, 
                  -908170874L, -652123112L, -221576806L, 1893421092L, -1337017914L, 
                  1810472933L, 1129781333L, 1414512861L, -1383904912L, 1553945052L, 
                  1959714339L, -737463321L, -1208684064L, -1806682705L, 274798480L, 
                  283185078L, 1394174670L, 1056174018L, 476966481L, 992957759L, 
                  1314393149L, 1626051874L, -1959865212L, 1934035236L, -1578656572L, 
                  1388860634L, -385346117L, -797213886L, -265885879L, -1369523083L, 
                  -109535459L, 63446241L, 142901344L, -869233916L, 1649283766L, 
                  -156001219L, -972161338L, -1785475312L, -2111545868L, -672239067L, 
                  -1793021768L, -1829083699L, 212149810L, -408582437L, 4398326L, 
                  -815423659L, 928236454L, 1558293485L, -1294777372L, 1687563396L, 
                  -1112158858L, 401013199L, 1118621393L, 1457662597L, 1758915456L, 
                  -1972203477L, 1616310339L, 568364066L, -2142105725L, 1113152765L, 
                  -880727309L, 463725441L, 68464333L, -1197337398L, 1535690804L, 
                  -336351145L, -1449160593L, -1777629838L, -193380907L, 128629040L, 
                  327722364L, 596765589L, -943329089L, 1178566603L, -1871500841L, 
                  -1693343807L, -648319107L, -4347453L, 258884252L, -1520999381L, 
                  -1328082201L, 1783130926L, 1423617136L, 940552L, 1854841677L, 
                  -808321327L, 1280062860L, 1712918026L, -1013656639L, 807777202L, 
                  -1454957880L, -1823314840L, 1701357697L, 539651538L, 193736663L, 
                  1146364944L, -1291420497L, -1651535721L, -1799524634L, 854016906L, 
                  -53837673L, -125343026L, 1174705402L, -621750140L, -953554677L, 
                  -1079572619L, 497705036L, -1834145204L, 1667496812L, 972695281L, 
                  230385907L, -496694749L, -1199420939L, -11293146L, 323874436L, 
                  -1462133397L, 879957117L, 2095594610L, -1671142743L)



data$w <- rep(0.5, N_obs_days)
data$w[data$cp_idx] <- 1

data$w[data$w == 0.5] <- 0.15

data$w[data$cp_idx] <- 0.2
data$w[c(88, 116, 142, 171, 175, 196, 229, 244, 247, 285, 293, 328, 358, 361)] <- 1

horseshoe_fit <- stan(file='stan_programs/fit_sales_horseshoe.stan',
                      data=data, seed=4938483, refresh=1000)

# Check diagnostics
capture.output(util$check_n_eff(horseshoe_fit))[1:5]
capture.output(util$check_rhat(horseshoe_fit))[1:5]
util$check_div(horseshoe_fit)
util$check_treedepth(horseshoe_fit)
util$check_energy(horseshoe_fit)

partition <- util$partition_div(horseshoe_fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

# Inner Core
par(mfrow=c(3, 3))

for (k in which(data$w == 0.5)[4*(1:9)]) {
  name_x <- paste("delta_days_tilde[", k, "]", sep='')
  name_y <- paste("lambda[", k, "]", sep='')
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("Day ", k),
       xlab=name_x, xlim=c(-1, 1), 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-8, 5))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
}

data$w[data$w == 0.5] <- 0.2

# Outer Core
par(mfrow=c(5, 5))

for (k in which(data$w == 1)) {
  name_x <- paste("delta_days_tilde[", k, "]", sep='')
  name_y <- paste("lambda[", k, "]", sep='')
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("Day ", k),
       xlab=name_x, xlim=c(-0.75, 0.75), 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-6, 7))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
}

data$w[data$cp_idx] <- 0.35
data$w[c(88, 116, 142, 171, 175, 196, 229, 244, 247, 285, 293, 328, 358, 361)] <- 1

# Try again
horseshoe_fit <- stan(file='stan_programs/fit_sales_horseshoe.stan',
                      data=data, seed=4938483, refresh=1000)

# Check diagnostics
util$check_all_diagnostics(horseshoe_fit)

partition <- util$partition_div(horseshoe_fit)
div_params <- partition[[1]]
nondiv_params <- partition[[2]]

# Inner Core
par(mfrow=c(3, 3))

for (k in which(data$w == 0.2)[4*(1:9)]) {
  name_x <- paste("delta_days_tilde[", k, "]", sep='')
  name_y <- paste("lambda[", k, "]", sep='')
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("Day ", k),
       xlab=name_x, xlim=c(-1, 1), 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-8, 5))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
}

# Outer Core
par(mfrow=c(5, 5))

for (k in which(data$w >= 0.35)) {
  name_x <- paste("delta_days_tilde[", k, "]", sep='')
  name_y <- paste("lambda[", k, "]", sep='')
  
  plot(nondiv_params[name_x][,1], log(nondiv_params[name_y][,1]),
       col=c_dark_trans, pch=16, main=paste("Day ", k),
       xlab=name_x, xlim=c(-0.75, 0.75), 
       ylab=paste("log(", name_y, ")", sep=""), ylim=c(-6, 7))
  points(div_params[name_x][,1], log(div_params[name_y][,1]),
         col=c_green_trans, pch=16)
}


data$w <- c(0.15, 0.15, 0.15, 0.15, 0.2, 0.15, 0.2, 0.15, 0.15, 0.15, 0.2, 
            0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.2, 0.15, 
            0.15, 0.2, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.15, 0.15, 0.2, 0.15, 0.15, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.2, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.2, 0.15, 0.2, 0.15, 0.15, 0.15, 0.2, 0.2, 
            0.15, 0.2, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 1, 
            0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.2, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.15, 0.15, 1, 0.15, 0.15, 0.15, 0.15, 0.15, 
            0.2, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 1, 0.15, 
            0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.2, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 
            0.2, 0.15, 0.15, 0.15, 0.15, 1, 0.15, 0.15, 0.15, 1, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.2, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.15, 0.2, 0.15, 0.2, 1, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.2, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.2, 0.15, 
            0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 1, 0.15, 0.15, 0.15, 
            0.15, 0.2, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.2, 0.15, 0.15, 
            1, 0.15, 0.15, 1, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 
            0.15, 0.2, 0.15, 0.15, 0.15, 0.15, 0.15, 1, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.15, 1, 0.15, 0.15, 0.15, 0.2, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.15, 0.2, 0.15, 0.2, 0.2, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.2, 0.2, 0.15, 0.2, 0.2, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 1, 0.15, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.2, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.2, 0.15, 0.15, 0.15, 0.2, 0.15, 0.2, 0.15, 
            0.15, 0.15, 0.15, 1, 0.15, 0.15, 1, 0.15, 0.15, 0.15, 0.15, 0.15, 
            0.15, 0.15, 0.15, 0.15)

.Random.seed <- c(10403L, 54L, -809302514L, -1649027925L, 836617953L, -1249118417L, 
                  1278836891L, 1685176797L, -1654588936L, 1905218982L, 1418063728L, 
                  -900832307L, 10935867L, 223600251L, 918001093L, -193246870L, 
                  -715013192L, -1147444885L, -48170971L, 287154586L, 898872228L, 
                  1985097474L, 813880538L, -1462121185L, 464361994L, -1679288726L, 
                  -1993959960L, 1986444395L, -1645791546L, 1218736261L, 104697543L, 
                  487323447L, 524642538L, -1186938912L, -1459728208L, 1573319732L, 
                  -922967540L, 496703551L, -2005034741L, -1162069396L, -588094218L, 
                  2084118371L, 1336334429L, 42159845L, -367497111L, -758974312L, 
                  -123439497L, -1445243774L, -1759287408L, -2090337883L, -1207465398L, 
                  512891L, -778291492L, 1508498582L, -1911834872L, 335919586L, 
                  -601280334L, -519957804L, 206864330L, -403042822L, 1769430502L, 
                  966491242L, -2142052785L, 639379867L, 964716993L, 1450175052L, 
                  -365018913L, -1209520458L, 1126038731L, -1197886125L, 1722141607L, 
                  -1199047730L, 1681976421L, -1876482649L, -425952786L, -834685866L, 
                  1730983462L, -1374057795L, -1300668353L, 2095553478L, -409587153L, 
                  1041759712L, -1505991315L, -1809682093L, 393386179L, 286872881L, 
                  1062790857L, -482369449L, -1490862928L, 1599068955L, 2077707579L, 
                  -218063015L, 556520195L, 117761664L, -815357251L, -153321107L, 
                  1395457942L, -418515380L, 1741531716L, 879268998L, 618519779L, 
                  665619798L, 952670250L, 297368708L, 1441631573L, -312543908L, 
                  -444007866L, -2084615288L, -1594160057L, -208111664L, 943488324L, 
                  1643677429L, -1442519946L, 79931351L, 1789585057L, 1008119529L, 
                  1477522689L, 959184526L, 1484627404L, -655905929L, 766565263L, 
                  -1800491409L, -1436807012L, 1902683776L, -1161109379L, 376785194L, 
                  1865022879L, 1645315816L, 1719529771L, 1302645265L, 2102266452L, 
                  -1965434073L, -120290649L, -118956406L, 1267690685L, -1001026116L, 
                  159582074L, 911065087L, 1679083262L, 1567882608L, 1848371673L, 
                  -1016171307L, 718969749L, 1134295365L, 271446877L, 2009938777L, 
                  -92941395L, -556135669L, -1849897848L, 1617656258L, 2053076467L, 
                  -1435652697L, 681414953L, 1993110684L, 1970802862L, -902581651L, 
                  -914370439L, -2122194650L, 1181431448L, -862623383L, -500929858L, 
                  -768624268L, 614616720L, 262822418L, 55117829L, 793387909L, -471556996L, 
                  -1617674827L, 926192080L, 2131335012L, -13819516L, -103189180L, 
                  -421002356L, 1168301701L, -1751731071L, -1573154676L, 1813326885L, 
                  -387334501L, 1112159171L, -220128778L, 1975424258L, 1264689871L, 
                  -997124153L, -2089682654L, -1623245438L, -1290180552L, -1221412066L, 
                  -778777468L, 2138238340L, -646798295L, -1239262575L, 552796697L, 
                  1348616574L, -144987778L, 1258715279L, -1936520233L, 2045187265L, 
                  2049723351L, -162145452L, 1479880025L, 1384696419L, 191385880L, 
                  1865399051L, -956415030L, 1582508754L, 1918936458L, -768794464L, 
                  -96937557L, 1396363941L, -394113433L, 765288709L, 1583549289L, 
                  2093964499L, 1643682059L, 2030422695L, -561671778L, 1806106102L, 
                  -636217839L, 778744373L, -1303863103L, 159629913L, 815279723L, 
                  -1324586728L, 1022926588L, -1999268082L, 1181263850L, -649907656L, 
                  -1256689846L, 1464015558L, -1184062827L, -223744786L, 2106675482L, 
                  858663020L, -1032182637L, -1808803529L, 1311778326L, -1139278779L, 
                  1721839891L, 1157423655L, 1000192085L, -1836182789L, -206233435L, 
                  -180870101L, 1716309376L, 989031072L, 1790442471L, 639569958L, 
                  384038419L, -977917626L, -792580150L, 1841803693L, -173293031L, 
                  -411274345L, -1776500492L, 1377798778L, -1077106283L, -349209973L, 
                  -774490401L, -1233572039L, -683888514L, -59702360L, 1407445834L, 
                  -949296052L, -980409565L, -1734407977L, 1141361783L, -41617754L, 
                  -1628720697L, 1403626520L, -889591515L, -262662905L, -1143335352L, 
                  116111341L, 701702832L, 2047137284L, 623236480L, 1093741154L, 
                  392118927L, -74217547L, 1025636891L, -142625745L, -1099822570L, 
                  1548730617L, 1965207878L, -216175698L, -728665762L, 517539541L, 
                  35751103L, 243980891L, -705264552L, -939374650L, -275884252L, 
                  939871047L, -495109493L, -1906814392L, 1323393329L, -105022075L, 
                  -1365165508L, -2056433065L, 467307105L, -2051658777L, 470078325L, 
                  -1133749909L, 2040299753L, 1177484256L, 1841931556L, -1676819337L, 
                  -911451489L, -1157751754L, 835177681L, 1643897980L, 1527957938L, 
                  2348978L, -881387386L, -1507999734L, -612486017L, 46847612L, 
                  -1799422625L, 2024118669L, -312907091L, -726467362L, 203820618L, 
                  1026147730L, 1989244774L, 1671372849L, 1878891750L, 1196138206L, 
                  622469620L, -1267854396L, 1392467275L, -56696252L, -853923739L, 
                  -2000398591L, -740453798L, 1246112720L, -794673496L, 1442680937L, 
                  783655076L, 700013672L, -1505782159L, 994550737L, -418180376L, 
                  -1966920885L, -365778063L, 488365994L, -556628251L, -2127004524L, 
                  -386401573L, 1745451968L, 1205113358L, 432773061L, -2117839630L, 
                  -557726508L, 1482643433L, -1651156357L, 1847675515L, 1316770L, 
                  1790574982L, 1854834157L, -1473940789L, 1794678962L, -355817063L, 
                  1217133699L, 1891774341L, -17942259L, 61987390L, -2129214625L, 
                  -846810315L, 376893477L, -1073314196L, 1262161900L, 717193001L, 
                  1100479742L, -1829011200L, -139601754L, -77521541L, -501967070L, 
                  1257155013L, -802734293L, -2098270454L, 827861875L, 1726788307L, 
                  -281554663L, 49888470L, 637741180L, -499082275L, -1217164996L, 
                  668990235L, -815615977L, 1945850457L, 1741817277L, 1851548601L, 
                  870124805L, -1439787676L, 1166421734L, -868138993L, 921203235L, 
                  1210522568L, -75112018L, 1304013213L, -872295368L, 1077866087L, 
                  1682260573L, -2108095741L, 520235509L, 1162003169L, -225664611L, 
                  -1372071407L, -2112165590L, -1979484341L, -1630877256L, -636848109L, 
                  1431591212L, -2123795028L, 843974784L, -532877801L, 953393990L, 
                  819740317L, 708479931L, 467554348L, -1227368880L, -999625313L, 
                  1605458723L, -881784612L, -1157183709L, 87900280L, 1573180896L, 
                  -1998122005L, 1755844923L, -293806485L, -1098544874L, -137597504L, 
                  -1222868541L, -590059898L, -2098292261L, -1524293286L, -120096091L, 
                  651930973L, 1550266628L, -29876983L, -1914481919L, -1314782060L, 
                  -2079006477L, 1104688370L, -1113441369L, -1628981824L, -1031558852L, 
                  424874107L, -1292723418L, 2135840786L, 58420209L, 534241037L, 
                  -841478802L, 429792318L, 203584717L, 1176050954L, 1669982184L, 
                  -1601595454L, 285923460L, 1164713879L, -699533839L, -1814358291L, 
                  -904191595L, -1649752477L, 1042742332L, 498534036L, -478231205L, 
                  -873795188L, -1427797050L, 1402337481L, 1702891913L, -1731273055L, 
                  762039641L, -933009058L, -889006415L, -1487199978L, -1769608932L, 
                  -1964574055L, -1457647348L, -82195063L, 342047643L, 127328441L, 
                  -1111107843L, -1316575983L, 1935500293L, 693429999L, 813124143L, 
                  -1470564320L, 1688073866L, 1894566102L, -759839983L, 892945585L, 
                  1743117242L, 226105754L, -1131174359L, -1641596741L, 1616028097L, 
                  143536006L, -360106430L, -219779721L, 1878453250L, 881461792L, 
                  1792173365L, -1460229752L, -1256843033L, 181839880L, -1714444504L, 
                  -908170874L, -652123112L, -221576806L, 1893421092L, -1337017914L, 
                  1810472933L, 1129781333L, 1414512861L, -1383904912L, 1553945052L, 
                  1959714339L, -737463321L, -1208684064L, -1806682705L, 274798480L, 
                  283185078L, 1394174670L, 1056174018L, 476966481L, 992957759L, 
                  1314393149L, 1626051874L, -1959865212L, 1934035236L, -1578656572L, 
                  1388860634L, -385346117L, -797213886L, -265885879L, -1369523083L, 
                  -109535459L, 63446241L, 142901344L, -869233916L, 1649283766L, 
                  -156001219L, -972161338L, -1785475312L, -2111545868L, -672239067L, 
                  -1793021768L, -1829083699L, 212149810L, -408582437L, 4398326L, 
                  -815423659L, 928236454L, 1558293485L, -1294777372L, 1687563396L, 
                  -1112158858L, 401013199L, 1118621393L, 1457662597L, 1758915456L, 
                  -1972203477L, 1616310339L, 568364066L, -2142105725L, 1113152765L, 
                  -880727309L, 463725441L, 68464333L, -1197337398L, 1535690804L, 
                  -336351145L, -1449160593L, -1777629838L, -193380907L, 128629040L, 
                  327722364L, 596765589L, -943329089L, 1178566603L, -1871500841L, 
                  -1693343807L, -648319107L, -4347453L, 258884252L, -1520999381L, 
                  -1328082201L, 1783130926L, 1423617136L, 940552L, 1854841677L, 
                  -808321327L, 1280062860L, 1712918026L, -1013656639L, 807777202L, 
                  -1454957880L, -1823314840L, 1701357697L, 539651538L, 193736663L, 
                  1146364944L, -1291420497L, -1651535721L, -1799524634L, 854016906L, 
                  -53837673L, -125343026L, 1174705402L, -621750140L, -953554677L, 
                  -1079572619L, 497705036L, -1834145204L, 1667496812L, 972695281L, 
                  230385907L, -496694749L, -1199420939L, -11293146L, 323874436L, 
                  -1462133397L, 879957117L, 2095594610L, -1671142743L)

# Charge!
horseshoe_fit <- stan(file='stan_programs/fit_sales_horseshoe.stan',
                      data=data, seed=4938483, refresh=1000, 
                      control=list(adapt_delta=0.999))

# Check diagnostics
util$check_all_diagnostics(horseshoe_fit)

# Ugh, close enough

# Check marginal posteriors
samples = extract(horseshoe_fit)

par(mfrow=c(2, 3))

hist(samples$alpha, breaks=seq(5.40, 5.45, 0.05 / 50),
     main="", xlab="alpha", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$alpha, col=c_light, lty=1, lw=2)

hist(samples$delta_season, seq(0.2, 0.3, 0.1 / 50),
     main="", xlab="delta_season", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$delta_season, col=c_light, lty=1, lw=2)

hist(samples$phi, seq(9.8, 10.2, 0.4 / 50), 
     main="", xlab="phi", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=truth$phi, col=c_light, lty=1, lw=2)

hist(samples$tau_day, seq(0, 0.1, 0.1 / 50),
     main="", xlab="tau_day", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)

hist(samples$milli_inv_psi, seq(0, 3, 3 / 50),
     main="", xlab="milli_inv_psi", yaxt='n', ylab="",
     col=c_dark, border=c_dark_highlight)
abline(v=1000 * truth$inv_psi, col=c_light, lty=1, lw=2)

par(mfrow=c(2, 1))

util$plot_excess_variation(samples, data, truth, "Horseshoe Population Model")

util$plot_excess_variation_residual(samples, data, truth, "Horseshoe Population Model")

infer_crps_comp["Horseshoe"] = util$sum_empirical_crps(samples$delta_days, truth$delta_days)
print(infer_crps_comp)

# Check expected sales time series
par(mfrow=c(2, 1))

util$plot_expected_sales(samples, data, "Horseshoe Population Model")

util$plot_expected_sales_residual(samples, data, "Horseshoe Population Model")

# Check posterior retrodictions and predictions
util$plot_dictions(samples, data, "Horseshoe Population Model")

util$plot_dictions_residual(samples, data, "Horseshoe Population Model")

# Evaluate predictive performance with cumulative risk probability score
pred_crps_comp["Horseshoe"] = util$sum_empirical_rps(samples$y_post_pre, data$y_pred)
print(pred_crps_comp)

############################################################
# Comparison
############################################################ 

par(mfrow=c(1, 1))

plot(as.numeric(infer_crps_comp[1,]), as.numeric(pred_crps_comp[1,]), 
     col=c_dark, pch=16, cex=1.5,
     xlim=c(0, 10), xlab="Inferential ECRPS",
     ylim=c(1000, 2500), ylab="Predictive ECRPS")

for (m in 1:5) {
  text(infer_crps_comp[1,m], pred_crps_comp[1,m], 
       cex=1.25, label=names(infer_crps_comp)[m], 
       pos=4, col=c_dark)
}
