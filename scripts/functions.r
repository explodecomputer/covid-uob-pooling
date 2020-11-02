library(tidyverse)
library(BBmisc)

simulate_viralload <- function(ids,f0,f1,f2,g0,g1,g2) {
  # lognormal distribution of viral load over time
  x <- seq(0.01, f0, by=0.01)
  y <- dnorm(log(x), mean=f1, sd=f2)
  a <- trapz(x,y)
  
  # scale maximum viral load
  l <- rbeta(nrow(ids),g0,g1)*g2
  ids$cvl <- a * (10^l)
  
  # sample each individual at a random time
  ts <- sample(1:length(x),nrow(ids),replace=TRUE)
  # ids$day <- x[ts]
  ids$vl <- y[ts] * (10^l)
  
  return(ids)
}

simulate_infection <- function(ids, prevalence, spread, containment)
{
	# Sample initial individuals
	dat <- tibble(
		id = sample(ids$id, nrow(ids) * prevalence, replace=FALSE),
		spread = rpois(length(id), spread)
		) %>%
		inner_join(ids, ., by="id")
	
	new_infecteds <- lapply(1:nrow(dat), function(i)
	{
		loc <- sample(1:3, dat$spread[i], replace=TRUE, prob=containment)
		return(c(
			sample(subset(ids, circle == dat$circle[i])$id, min(dat$count[i], sum(loc==1)), replace=FALSE),
			sample(subset(ids, location == dat$location[i])$id, sum(loc==2), replace=FALSE),
			sample(ids$id, sum(loc==3), replace=FALSE)
		))
	}) %>% unlist %>% unique
	ids$infected <- ids$id %in% c(dat$id, new_infecteds)
	ids <- ids %>%
		group_by(location) %>%
		mutate(location_outbreak = sum(infected >= 2)) %>%
		ungroup()
	
	# remove vl for anyone not infected
	ids$vl[ids$infected==FALSE] <- 0
	# ids$cvl[ids$infected==FALSE] <- 0
	
	return(ids)
}

allocate_pools <- function(ids, pool_size, random=FALSE)
{
	random_pools <- function(x, pool_size)
	{
		npool <- ceiling(length(x) / pool_size)
		rep(1:npool, each=pool_size)[1:length(x)] %>% sample
	}

	if(random)
	{
		ids$assay_pool <- random_pools(1:nrow(ids), pool_size)
		ids <- ids %>%
			group_by(assay_pool) %>%
			mutate(poolsize=n()) %>%
			ungroup()
		return(ids)
	}

	# First break up large circles 
	temp1 <- subset(ids, count > pool_size) %>%
		group_by(circle) %>%
		mutate(
			poolbreak = cut(1:n(), ceiling(n()/pool_size)) %>% as.numeric(),
			pool = paste0(circle, "-", poolbreak)
		) %>% ungroup()
	temp2 <- subset(ids, count <= pool_size) %>%
		mutate(
			poolbreak = 1,
			pool = paste0(circle, "-", poolbreak)
		)
	ids <- bind_rows(temp1, temp2)
	dat <- ids %>%
		group_by(pool) %>%
		summarise(poolsize=n(), .groups="drop")
	# Note that `cut` tries to make even cuts. this might not be optimal
	# Use bin packing algorithm to assign to pools
	dat$assay_pool <- BBmisc::binPack(dat$poolsize, pool_size)
	ids <- inner_join(ids, dat, by="pool") %>%
		group_by(assay_pool) %>%
		mutate(poolsize = n()) %>% ungroup()
	return(ids)
}

assay_sensitivity_quad <- function(dilution, baseline, beta)
{
	vec <- -beta * (dilution-1)^2 + baseline
	vec[vec < 0] <- 0
	vec
}

assay_sensitivity_decay <- function(dilution, baseline, beta)
{
	exp(-(dilution-1) * beta) * baseline
}

assay_sensitivity_sigmoid <- function(dilution)
{
	1 / (1 + exp(-dilution/4+4)) * -1 + 1
}


simulate_testing <- function(ids)
{
	detected_pools <- ids %>% group_by(assay_pool) %>%
		summarise(
			ninfected = sum(infected), 
			sensitivity = assay_sensitivity_sigmoid(n()/ninfected),
			pool_detected = rbinom(1, 1, sensitivity),
			.groups="drop"
		)
	baseline_sensitivity <- assay_sensitivity_sigmoid(1)
	ids <- inner_join(ids, detected_pools, by="assay_pool")
	ids <- ids %>% group_by(circle) %>%
		mutate(circle_detected = any(pool_detected))
	ids$id_detected <- rbinom(nrow(ids), 1, baseline_sensitivity * ids$infected)
	ids <- ids %>% 
		group_by(location) %>%
		mutate(
			ind_outbreak_detected = sum(id_detected) > 2,
			pool_outbreak_detected = sum(pool_detected) > 2,
			poolfollowup_outbreak_detected = sum(pool_detected * id_detected) > 2
		) %>% ungroup()
	return(ids)
}


summarise_simulations <- function(ids, cost_samplingkit, cost_test)
{
	npool <- length(unique(ids$assay_pool))
	nfollowup <- sum(ids$pool_detected)

	ob <- group_by(ids, location) %>%
		summarise(
			outbreak = any(location_outbreak), 
			ind_outbreak_detected = any(ind_outbreak_detected),
			pool_outbreak_detected = any(pool_outbreak_detected),
			poolfollowup_outbreak_detected = any(poolfollowup_outbreak_detected),
			.groups="drop"
		)

	tibble(
		nstudents = nrow(ids),
		true_prevalence = sum(ids$infected) / nrow(ids),
		prevalence.ind_tests = sum(ids$id_detected) / nstudents,
		prevalence.pool_tests = sum(ids$pool_detected) / nstudents,
		prevalence.poolfollowup_tests = sum(ids$pool_detected & ids$id_detected) / nstudents,

		cost.ind_tests = nrow(ids) * cost_samplingkit + nrow(ids) * cost_test,
		cost.pool_tests = nrow(ids) * cost_samplingkit + npool * cost_test,
		cost.poolfollowup_tests = cost.pool_tests + nfollowup * cost_test,

		reagentuse.ind_tests = nrow(ids),
		reagentuse.pool_tests = npool,
		reagentuse.poolfollowup_tests = npool + nfollowup,

		sensitivity.ind_tests = sum(ids$id_detected) / sum(ids$infected),
		sensitivity.pool_tests = sum(ids$pool_detected * ids$infected) / sum(ids$infected),
		sensitivity.poolfollowup_tests = sum(ids$pool_detected * ids$infected * ids$id_detected) / sum(ids$infected),

		specificity.ind_tests = sum(!ids$id_detected & !ids$infected) / sum(!ids$infected),
		specificity.pool_tests = sum(!ids$pool_detected & !ids$infected) / sum(!ids$infected),
		specificity.poolfollowup_tests = 1 - sum(ids$pool_detected & ids$id_detected & !ids$infected) / sum(!ids$infected),

		sensitivity.outbreak.ind_tests = sum(ob$ind_outbreak_detected) / sum(ob$outbreak),
		sensitivity.outbreak.pool_outbreak_detected = sum(ob$pool_outbreak_detected) / sum(ob$outbreak),
		sensitivity.outbreak.poolfollowup_outbreak_detected = sum(ob$poolfollowup_outbreak_detected) / sum(ob$outbreak),

		specificity.outbreak.ind_tests = sum(!ob$ind_outbreak_detected & !ob$outbreak) / sum(!ob$outbreak),
		specificity.outbreak.pool_outbreak_detected = sum(!ob$pool_outbreak_detected & !ob$outbreak) / sum(!ob$outbreak),
		specificity.outbreak.poolfollowup_outbreak_detected = sum(!ob$poolfollowup_outbreak_detected & !ob$outbreak) / sum(!ob$outbreak)
	)
}


run_simulation <- function(ids, param)
{
	stopifnot(nrow(param) == 1)
	if(param$containment == "high")
	{
		containment <- c(0.9, 0.09, 0.01)
	} else if(param$containment == "medium") {
		containment <- c(0.6, 0.3, 0.1)
	} else if(param$containment == "low") {
		containment <- c(0.34, 0.33, 0.33)
	} else {
		stop("containment value")
	}
  
  ids <- simulate_viralload(ids,param$f0,param$f1,param$f2,param$g0,param$g1,param$g2) ## 
	ids <- simulate_infection(ids, param$prevalence, spread=param$spread, containment)
	ids <- allocate_pools(ids, param$pool_size, random=param$random_pooling)
	ids <- simulate_testing(ids)
	res <- summarise_simulations(ids, param$cost_samplingkit, param$cost_test)
	return(bind_cols(param, res))
}



