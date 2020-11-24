library(parallel)
library(dplyr)

###

source("functions.r")
load("../data/circles.rdata")
load("../data/efficiency_params.rdata")

params <- expand.grid(
	# infection characteristics
	prevalence = c(0.001, 0.01, 0.05),
	spread = c(0.5, 1, 3),
	containment = c("high", "medium", "low"),

	# pooling characteristic
	pool_size = c(2, 3, 4, 5, 10, 15, 20, 25, 30),
	random_pooling = c(TRUE, FALSE),

	# costs
	cost_samplingkit = 2,
	cost_test = 25,
	replicates = c(1:100)
)

params$f0 <- 14 # number of days of viral load
params$f1 <- 0 # mean of lognormal for viral load distribution
params$f2 <- 0.8 # sd of lognormal for viral load distribution
params$g0 <- 0.3 # beta shape 1
params$g1 <- 1 # beta shape 2
params$g2 <- 3 # multiplier for beta distribution
params$Emin <- 0.1 # Minimum PCR efficiency
params$Emax <- 0.95 # Maximum PCR efficiency
params$Ea <- efficiency_params$par[1] # Beta distribution a parameter for efficiency
params$Eb <- efficiency_params$par[2] # Beta distribution b parameter for efficiency
params$Ct <- 35 # Number cycles for detection
params$Rct <- Rct # Log Rct fluourescence detection value
params$fp <- 0.005 # Testing false positive rate (per test)

res <- mclapply(1:nrow(params), function(i) {
	message(i, " of ", nrow(params))
	run_simulation(ids, params[i,])
}, mc.cores=16) %>% bind_rows()

save(res, file="../results/sim.rdata")
