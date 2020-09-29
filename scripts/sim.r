library(parallel)
library(dplyr)

###

source("functions.r")
load("../data/circles.rdata")

params <- expand.grid(
	# infection characteristics
	prevalence = c(0.001, 0.01, 0.05),
	spread = c(0.5, 1, 2),
	containment = c("high", "medium", "low"),

	# pooling characteristic
	pool_size = c(2, 3, 4, 5, 10, 15, 20, 25, 30),
	random_pooling = c(TRUE, FALSE),

	# costs
	cost_samplingkit = 2,
	cost_test = 25,
	replicates = c(1:100)
)

res <- mclapply(1:nrow(params), function(i) {
	message(i, " of ", nrow(params))
	run_simulation(ids, params[i,])
}, mc.cores=16) %>% bind_rows()

save(res, file="../results/sim.rdata")
