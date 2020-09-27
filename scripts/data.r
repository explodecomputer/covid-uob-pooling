library(tidyverse)

dat1 <- readxl::read_xlsx("../data/Living Circles Count Update.xlsx", sheet=1, skip=2)
names(dat1) <- c("circle", "count", "revisedcount", "rawcount", "communitycircle")
dat2 <- readxl::read_xlsx("../data/Living Circles Count Update.xlsx", sheet=2)
names(dat2) <- c("location", "circle", "description", "type")
dat2 <- dat2 %>%
	subset(., !duplicated(circle))
circles <- dplyr::inner_join(dat1, dat2) %>%
	dplyr::select(location, circle, count)

# Create ids in each circle
ids <- lapply(1:nrow(circles), function(i)
	{
		tibble(
			location=circles$location[i],
			circle=circles$circle[i],
			id=paste0(circles$circle[i], ":", 1:circles$count[i]),
			count=circles$count[i]
		)
	}) %>% bind_rows()


save(circles, ids, file="../data/circles.rdata")
