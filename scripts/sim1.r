library(parallel)
library(dplyr)

###

source("functions.r")
load("../data/circles.rdata")

params <- expand.grid(
	# infection characteristics
	prevalence = c(0.001, 0.005, 0.01, 0.05),
	spread = c(0.5, 1, 2),
	containment = c("high", "medium", "low"),

	# pooling characteristic
	pool_size = c(3, 5, 8, 10, 15, 20, 25, 30),
	
	# assay characteristics
	baseline_sensitivity = 0.98,
	dilution_decay = c(0.05, 0.1, 0.2),

	# costs
	cost_samplingkit = 1,
	cost_extraction = 1,
	cost_pcr = 1,
	replicates = c(1:1)
)

res <- mclapply(1:nrow(params), function(i) {
	message(i, " of ", nrow(params))
	run_simulation(ids, params[i,])
}, mc.cores=16) %>% bind_rows()

save(res, file="../data/sim1.rdata")
