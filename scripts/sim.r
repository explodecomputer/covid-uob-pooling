library(parallel)
library(dplyr)

###

source("functions.r")
load("../data/circles.rdata")
load("../data/efficiency_params.rdata")
load("../data/lfd_fit.rdata")
load("../data/containment.rdata")


params <- expand.grid(
	# infection characteristics
	prevalence = c(0.001, 0.01, 0.05),
	spread = c(0.8, 1, 3),
	containment = c("high", "conquest", "scs1"),

	# pooling characteristic
	pool_size = c(2, 3, 4, 5, 10, 15, 20, 25, 30),
	random_pooling = c(TRUE, FALSE),

	# costs
	cost_samplingkit = 2,
	cost_test = 25,
	replicates = c(1:100)
)

params$ct_max <- efficiency_params$par[1] # maximum ct value
params$ct_min <- efficiency_params$par[2] # minimum ct value
params$ct_alpha <- efficiency_params$par[3] # Beta distribution a parameter for ct
params$ct_beta <- efficiency_params$par[4] # Beta distribution b parameter for ct
params$e_alpha <- efficiency_params$par[5] # Beta distribution a parameter for efficiency
params$e_beta <- efficiency_params$par[6] # Beta distribution b parameter for efficiency
params$e_min <- 0.65 # Minimum PCR efficiency
params$e_max <- 0.9 # Maximum PCR efficiency
params$ctthresh <- 35 # Number cycles for detection
params$rct <- 10 # Log Rct fluourescence detection value - arbitrary
params$pcr_fp <- 0.005 # Testing false positive rate (per test)
params$lfd_Asym <- lfd_fit$coef[1,1]
params$lfd_xmid <- lfd_fit$coef[2,1]
params$lfd_scal <- lfd_fit$coef[3,1]
params$lfd_fp <- 0.0032
params$lfd_cost <- 5

res <- mclapply(1:nrow(params), function(i) {
	message(i, " of ", nrow(params))
	run_simulation(ids, params[i,], containment)
}, mc.cores=16) %>% bind_rows()

save(res, file="../results/sim.rdata")
