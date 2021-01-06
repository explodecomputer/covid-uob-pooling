library(tidyverse)
library(BBmisc)
library(caTools)

simulate_viral_load <- function(ids, ct_alpha, ct_beta, ct_max, ct_min, e_alpha, e_beta, e_max, e_min, rct, gamma_shape=2, gamma_scale=2)
{
	n <- nrow(ids)
	ct <- rbeta(n, ct_alpha, ct_beta) * (ct_max - ct_min) + ct_min
	e <- rbeta(n, e_alpha, e_beta) * (e_max-e_min) + e_min
	vl <- rct / (1+e)^ct
	ids$id_E <- e
	ids$id_ct <- ct
	ids$id_vl <- vl

	vlshape <- tibble(
		t=seq(0,21,length.out=2101),
		vlshape=dgamma(t, shape=gamma_shape, scale=gamma_scale)
	)
	vlshape$vlshape <- vlshape$vlshape / max(vlshape$vlshape)
	len <- nrow(subset(vlshape, t <= 1))
	i1 <- sample(1:(which(vlshape$t == 18)[1]), n, replace=TRUE)
	i2 <- i1 + len * 3
	max_vl <- rct / (1 + e)^ct

	vl1 <- max_vl * vlshape$vlshape[i1]
	vl2 <- max_vl * vlshape$vlshape[i2]

	ct1 <- log(rct/vl1) / log(1+e)
	ct2 <- log(rct/vl2) / log(1+e)

	ids$id_ct1 <- ct1
	ids$id_ct2 <- ct2
	ids$id_vl1 <- vl1
	ids$id_vl2 <- vl2
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
	ids$id_vl[ids$infected==FALSE] <- 0
	ids$id_vl1[ids$infected==FALSE] <- 0
	ids$id_vl2[ids$infected==FALSE] <- 0
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

assay_ppv <- function(sensitivity, prevalence, specificity)
{
	sensitivity * prevalence / (sensitivity * prevalence + (1-specificity) * (1-prevalence))
}

lfd_prob <- function(x, asym, xmid, scal)
{
	asym / (1 + exp((xmid - x)/scal))
}

simulate_testing <- function(ids, ctthresh, rct, pcr_fp, lfd_Asym, lfd_xmid, lfd_scal, lfd_fp)
{
	ids$id_detected <- ids$id_ct < ctthresh
	ids$id_detected1 <- ids$id_ct1 < ctthresh
	ids$id_detected2 <- ids$id_ct2 < ctthresh
	ids$id_detected[!ids$infected] <- rbinom(sum(!ids$infected), 1, pcr_fp)
	ids$id_detected1[!ids$infected] <- rbinom(sum(!ids$infected), 1, pcr_fp)
	ids$id_detected2[!ids$infected] <- rbinom(sum(!ids$infected), 1, pcr_fp)
	ids$id_lfd_detected <- rbinom(nrow(ids), 1, lfd_prob(ids$id_ct, lfd_Asym, lfd_xmid, lfd_scal))
	ids$id_lfd_detected1 <- rbinom(nrow(ids), 1, lfd_prob(ids$id_ct1, lfd_Asym, lfd_xmid, lfd_scal))
	ids$id_lfd_detected2 <- rbinom(nrow(ids), 1, lfd_prob(ids$id_ct2, lfd_Asym, lfd_xmid, lfd_scal))
	ids$id_lfd_detected[!ids$infected] <- rbinom(sum(!ids$infected), 1, lfd_fp)
	ids$id_lfd_detected1[!ids$infected] <- rbinom(sum(!ids$infected), 1, lfd_fp)
	ids$id_lfd_detected2[!ids$infected] <- rbinom(sum(!ids$infected), 1, lfd_fp)

	detected_pools <- ids %>% group_by(assay_pool) %>%
		summarise(
			pool_ninfected = sum(infected), 
			pool_infected = pool_ninfected > 0,
			pool_E = min(id_E),
			pool_vl = sum(id_vl),
			pool_vl1 = sum(id_vl1),
			pool_vl2 = sum(id_vl2),
			pool_Rn = pool_vl / n() * (1 + pool_E)^ctthresh,
			pool_Rn1 = pool_vl1 / n() * (1 + pool_E)^ctthresh,
			pool_Rn2 = pool_vl2 / n() * (1 + pool_E)^ctthresh,
			pool_ct = log(rct / pool_vl) / log((1 + pool_E)),
			pool_ct1 = log(rct / pool_vl1) / log((1 + pool_E)),
			pool_ct2 = log(rct / pool_vl2) / log((1 + pool_E)),
			pool_detected = pool_ct < ctthresh,
			pool_detected1 = pool_ct1 < ctthresh,
			pool_detected2 = pool_ct2 < ctthresh,
			.groups="drop"
		)
	detected_pools$pool_detected[!detected_pools$pool_infected] <- rbinom(sum(!detected_pools$pool_infected), 1, pcr_fp)
	detected_pools$pool_detected1[!detected_pools$pool_infected] <- rbinom(sum(!detected_pools$pool_infected), 1, pcr_fp)
	detected_pools$pool_detected2[!detected_pools$pool_infected] <- rbinom(sum(!detected_pools$pool_infected), 1, pcr_fp)
	ids <- inner_join(ids, detected_pools, by="assay_pool")
	ids <- ids %>% group_by(circle) %>%
		mutate(
			circle_detected = any(pool_detected),
			circle_detected1 = any(pool_detected1),
			circle_detected2 = any(pool_detected2)
		)
	ids <- ids %>% 
		group_by(location) %>%
		mutate(
			ind_outbreak_detected = sum(id_detected) > 2,
			pool_outbreak_detected = sum(pool_detected) > 2,
			poolfollowup_outbreak_detected = sum(pool_detected * id_detected) > 2
		) %>% ungroup()
	return(ids)
}


summarise_simulations <- function(ids, cost_samplingkit, cost_test, cost_lfd)
{
	print(head(ids))
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
		cost.lfd = nrow(ids) * cost_lfd,

		reagentuse.ind_tests = nrow(ids),
		reagentuse.pool_tests = npool,
		reagentuse.poolfollowup_tests = npool + nfollowup,

		sensitivity.ind_tests = sum(ids$infected & ids$id_detected) / sum(ids$infected),
		sensitivity.ind_tests1 = sum(ids$infected & ids$id_detected1) / sum(ids$infected),
		sensitivity.ind_tests2 = sum(ids$infected & ids$id_detected2) / sum(ids$infected),
		sensitivity.pool_tests = sum(ids$pool_detected * ids$infected) / sum(ids$infected),
		sensitivity.pool_tests1 = sum(ids$pool_detected1 * ids$infected) / sum(ids$infected),
		sensitivity.pool_tests2 = sum(ids$pool_detected2 * ids$infected) / sum(ids$infected),
		sensitivity.poolfollowup_tests = sum(ids$pool_detected * ids$infected * ids$id_detected) / sum(ids$infected),
		sensitivity.poolfollowup_tests1 = sum(ids$pool_detected1 * ids$infected * ids$id_detected1) / sum(ids$infected),
		sensitivity.poolfollowup_tests2 = sum(ids$pool_detected2 * ids$infected * ids$id_detected2) / sum(ids$infected),
		sensitivity.lfd = sum(ids$infected & ids$id_lfd_detected) / sum(ids$infected),
		sensitivity.lfd1 = sum(ids$infected & ids$id_lfd_detected1) / sum(ids$infected),
		sensitivity.lfd2 = sum(ids$infected & ids$id_lfd_detected2) / sum(ids$infected),
		sensitivity.lfd_ever = sum(ids$infected & (ids$id_lfd_detected1 | ids$id_lfd_detected2)) / sum(ids$infected),

		specificity.ind_tests = sum(!ids$id_detected & !ids$infected) / sum(!ids$infected),
		specificity.ind_tests1 = sum(!ids$id_detected1 & !ids$infected) / sum(!ids$infected),
		specificity.ind_tests2 = sum(!ids$id_detected2 & !ids$infected) / sum(!ids$infected),
		specificity.pool_tests = sum(!ids$pool_detected & !ids$infected) / sum(!ids$infected),
		specificity.pool_tests1 = sum(!ids$pool_detected1 & !ids$infected) / sum(!ids$infected),
		specificity.pool_tests2 = sum(!ids$pool_detected2 & !ids$infected) / sum(!ids$infected),
		specificity.poolfollowup_tests = 1 - sum(ids$pool_detected & ids$id_detected & !ids$infected) / sum(!ids$infected),
		specificity.poolfollowup_tests1 = 1 - sum(ids$pool_detected1 & ids$id_detected1 & !ids$infected) / sum(!ids$infected),
		specificity.poolfollowup_tests2 = 1 - sum(ids$pool_detected2 & ids$id_detected2 & !ids$infected) / sum(!ids$infected),
		specificity.lfd = sum(!ids$id_lfd_detected & !ids$infected) / sum(!ids$infected),
		specificity.lfd1 = sum(!ids$id_lfd_detected1 & !ids$infected) / sum(!ids$infected),
		specificity.lfd2 = sum(!ids$id_lfd_detected2 & !ids$infected) / sum(!ids$infected),
		specificity.lfd_ever = sum(!(ids$id_lfd_detected1 | ids$id_lfd_detected2) & !ids$infected) / sum(!ids$infected),

		sensitivity.outbreak.ind_tests = sum(ob$ind_outbreak_detected) / sum(ob$outbreak),
		sensitivity.outbreak.pool_outbreak_detected = sum(ob$pool_outbreak_detected) / sum(ob$outbreak),
		sensitivity.outbreak.poolfollowup_outbreak_detected = sum(ob$poolfollowup_outbreak_detected) / sum(ob$outbreak),

		specificity.outbreak.ind_tests = sum(!ob$ind_outbreak_detected & !ob$outbreak) / sum(!ob$outbreak),
		specificity.outbreak.pool_outbreak_detected = sum(!ob$pool_outbreak_detected & !ob$outbreak) / sum(!ob$outbreak),
		specificity.outbreak.poolfollowup_outbreak_detected = sum(!ob$poolfollowup_outbreak_detected & !ob$outbreak) / sum(!ob$outbreak)
	)
}


run_simulation <- function(ids, param, containment)
{
	stopifnot(nrow(param) == 1)
	stopifnot(param$containment %in% names(containment))
	cont <- containment[[param$containment]]

	ids <- simulate_viral_load(ids, param$ct_alpha, param$ct_beta, param$ct_max, param$ct_min, param$e_alpha, param$e_beta, param$e_max, param$e_min, param$rct)
	ids <- simulate_infection(ids, param$prevalence, spread=param$spread, cont)
	ids <- allocate_pools(ids, param$pool_size, random=param$random_pooling)
	ids <- simulate_testing(ids, param$ctthresh, param$rct, param$pcr_fp, param$lfd_Asym, param$lfd_xmid, param$lfd_scal, param$lfd_fp)
	res <- summarise_simulations(ids, param$cost_samplingkit, param$cost_test, param$lfd_cost)
	return(bind_cols(param, res))
}
